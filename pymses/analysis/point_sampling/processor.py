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
"""

"""
from ..data_processor import DataProcessor
try:
    from ..data_processor import DataProcess
except ImportError:
    DataProcess = object

import numpy as N


def point_sampling_process_task(idomain, amr_source, points_datamap, points, sample_point_kwargs, verbose):
    # Select points from current domain only
    ipoint_domain = N.nonzero(points_datamap == idomain)[0]
    points_in_domain = points[ipoint_domain, :]

    # Get tha AMR dataset #idomain
    amr_dset = amr_source.get_domain_dset(idomain, verbose)

    # Sample AMR fields in the domain #idomain
    points_dset = amr_dset.sample_points(points_in_domain, **sample_point_kwargs)

    return ipoint_domain, points_dset


class PointSamplingProcess(DataProcess):
    def __init__(self, queues, amr_source, points_datamap, points, sample_point_kwargs, verbose):
        tasks_queue, results_queue = queues
        super(PointSamplingProcess, self).__init__(tasks_queue, results_queue)
        self.amr_source = amr_source
        self.points_datamap = points_datamap
        self.points = points
        self.sample_point_kwargs = sample_point_kwargs
        self.ipoint_domain_list = []
        self.points_dset_list = []
        self.verbose = verbose

    def process_task(self, idomain):
        task_res = point_sampling_process_task(idomain, self.amr_source, self.points_datamap, self.points,
                                               self.sample_point_kwargs, self.verbose)

        self.ipoint_domain_list.append(task_res[0])
        self.points_dset_list.append(task_res[1])
        return None

    @property
    def result(self):
        res = (self.ipoint_domain_list, self.points_dset_list)
        return res


class PointSamplingProcessor(DataProcessor):
    r"""
    Create point-based data from AMR-based data by point sampling. Samples all available fields
    of the `amr_source` at the coordinates of the `points`.

    Parameters
    ----------
    amr_source : :class:`~pymses.sources.ramses.output.RamsesAmrSource`
        AMR data source
    verbose: ``bool``
        verbosity boolean flag.
    """
    def __init__(self, amr_source, verbose=None):
        super(PointSamplingProcessor, self).__init__(amr_source, None, amr_mandatory=True, verbose=verbose)
        self._points_datamap = None
        self._points = None
        self._sample_points_kwargs = {}

    def process_factory(self, tasks_queue, results_queue, verbosity):
        return PointSamplingProcess((tasks_queue, results_queue), self._source, self._points_datamap, self._points,
                                     self._sample_points_kwargs, verbosity)

    def process(self, points, add_cell_center=False, add_level=False, max_search_level=None,
                interpolation=False):
        """
        Point-sampling data processor processing method.

        Parameters
        ----------
        points : (`npoints`, `ndim`) ``array``
            sampling points coordinates
        add_level : ``boolean`` (default False)
            whether we need to add a `level` field in the returned dataset containing
            the value of the AMR level the sampling points fall into
        add_cell_center : ``boolean`` (default False)
            whether we need to add a `cell_center` field in the returned dataset containing
            the coordinates of the AMR cell center the sampling points fall into
        interpolation : ``boolean`` (default False)
            Experimental : A proper bi/tri-linear interpolation could be great!
            THIS IS NOT IMPLEMENTED YET : in this attempt we supposed corner cell data
            while ramses use centered cell data, letting alone the problem
            of different AMR level...

        Returns
        -------
        dset : :class:`~pymses.core.datasets.PointDataset`
            Contains all these sampled values.
        """
        self._points = N.asarray(points)
        if max_search_level is not None and not isinstance(max_search_level, int):
            raise AttributeError("'max_search_level' must be an integer value or None (default).")

        # Compute the domain ids of the points
        if self._source.dom_decomp is not None:
            self._points_datamap = self._source.dom_decomp.map_points(self._points)
        else:
            self._points_datamap = N.ones(self._points.shape[0])
        # Sort the points by datamap
        unique_domains = N.unique(self._points_datamap)
        ndomains = unique_domains.size

        ipoint_batches = []
        point_dsets = []

        # Fill sample_points option dict.
        # Set the AMR data source max. level to read
        if max_search_level is not None:
            self._source.set_read_levelmax(max_search_level)

        self._sample_points_kwargs["add_cell_center"] = add_cell_center
        self._sample_points_kwargs["add_level"] = add_level
        self._sample_points_kwargs["max_search_level"] = max_search_level
        self._sample_points_kwargs["interpolation"] = interpolation

        if self.use_multiprocessing and ndomains >= 2:
            self.init_multiprocessing_run(ndomains)

            # Loop over unique domains : enqueue jobs
            for idomain in unique_domains:
                self.enqueue_task(idomain)

            for res in self.get_results():
                ipoint_batches += res[0]
                point_dsets += res[1]
        else:
            # Loop over unique domains
            for idomain in unique_domains:
                task_res = point_sampling_process_task(idomain, self._source, self._points_datamap, self._points,
                                                       self._sample_points_kwargs, self._verbose)

                ipoint_batches.append(task_res[0])
                point_dsets.append(task_res[1])

        # Perform final reordering concatenation
        cat_method = point_dsets[0].concatenate
        return cat_method(point_dsets, reorder_indices=ipoint_batches)


__all__ = ["PointSamplingProcessor"]
