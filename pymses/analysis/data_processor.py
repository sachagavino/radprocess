# -*- coding: utf-8 -*-
#   This file is part of PyMSES.
#
#   PyMSES is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyMSES is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.

from pymses import rcConfig as pymsesrc
from pymses.core import DataSource

try:
    from multiprocessing import Process, Pipe, Queue, JoinableQueue, cpu_count
    _multiproc__mod_available = True
except ImportError:
    Process = object
    _multiproc__mod_available = False


class DataProcess(Process):
    poison_pill = None

    def __init__(self, task_queue, result_queue):
        super(DataProcess, self).__init__()
        self._task_queue = task_queue
        self._result_queue = result_queue

    def process_task(self, task):
        raise NotImplementedError()

    @property
    def result(self):
        raise NotImplementedError()

    def run(self):
        proc_name = self.name
        while True:
            next = self._task_queue.get()
            if next is DataProcess.poison_pill:
                # print "Poison pill #%s" % proc_name
                self._result_queue.put(self.result)
                # Poison pill induces process shutdown
                self._task_queue.task_done()
                break

            ret = self.process_task(next)
            # print "Processed task #%d" % next
            if ret is not None:
                # Reenqueue following task
                self._task_queue.put(ret)
            self._task_queue.task_done()


class DataProcessor(object):
    r"""
    Data processing abstract class

    Parameters
    ----------
    source : :class:`~pymses.core.sources.DataSource`
        data source
    operator: :class:`~pymses.analysis.operator.AbstractOperator`
        physical quantity data operator
    amr_mandatory: ``bool``
        Is an AMR data source mandatory for this DataProcessor instance. Default False.
    verbose: ``bool``
        verbosity boolean flag. Default None
    """
    def __init__(self, source, operator, amr_mandatory=False, verbose=None):
        if amr_mandatory and source.source_type() != DataSource.AMR_SOURCE:
            raise AttributeError("Data source must be of type 'AMR' !")
        self._source = source
        self._operator = operator
        self._user_needs_multiprocessing = True
        self._processes = None
        self._tasks_q = None
        self._results_q = None
        self._verbose = verbose

    def process(self, *args, **kwargs):
        r"""
        Processing abstract method

        Parameters
        ----------
        args : ``list``
            list of positional arguments
        kwargs : ``dict``
            keyword argument dictionary
        """
        raise NotImplementedError()

    def disable_multiprocessing(self):
        self._user_needs_multiprocessing = False

    @property
    def use_multiprocessing(self):
        if not _multiproc__mod_available:
            return False

        if cpu_count() == 1 or pymsesrc.multiprocessing_max_nproc <= 1:
            return False

        return self._user_needs_multiprocessing

    @classmethod
    def ncpu_pool(cls, required=None):
        if not _multiproc__mod_available:
            return 1

        min_avail = min(cpu_count(), pymsesrc.multiprocessing_max_nproc)
        if required is not None:
            return min(required, min_avail)
        else:
            return min_avail

    def process_factory(self, tasks_queue, results_queue, verbosity):
        raise NotImplementedError()

    def init_multiprocessing_run(self, ntasks):
        """
        TODO

        Parameters
        ----------
        ntasks: ``int``

        Returns
        -------

        """
        if not self.use_multiprocessing:
            return

        # Cleanup
        if self._processes is not None:
            for iproc in range(len(self._processes)-1, -1, -1):
                del self._processes[iproc]
        del self._processes

        # Establish communication queues
        self._tasks_q = JoinableQueue()
        self._results_q = Queue()

        # Start consumers
        self._nproc = self.ncpu_pool(required=ntasks)
        self._processes = [self.process_factory(self._tasks_q, self._results_q, self._verbose)
                           for i in range(self._nproc)]

        for dp in self._processes:
            dp.start()

    def enqueue_task(self, task):
        """
        Add a task object/index to the task queue

        task: task object or index to add to task joinable queue
        """
        self._tasks_q.put(task)

    def get_results(self):
        """
        TODO

        Returns
        -------
        """
        add_poison_pills = True
        if add_poison_pills:
            # Add a poison pill for each process
            for i in range(self._nproc):
                self._tasks_q.put(DataProcess.poison_pill)

        # Wait for all of the tasks to finish
        self._tasks_q.join()

        # Start yielding results
        for i in range(self._nproc):
            yield self._results_q.get()


__all__ = ["DataProcessor", 'DataProcess']
