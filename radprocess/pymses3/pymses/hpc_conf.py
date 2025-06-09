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
High-performance computing configuration package (MPI framework setup)

"""
try:
    from mpi4py import MPI as mpi4py_MPI
    mpi4py_mpi_comm_world = mpi4py_MPI.COMM_WORLD
except ImportError:
    mpi4py_MPI = None
    mpi4py_mpi_comm_world = None


class HighPerformanceComputingConfiguration(object):
    OPERATION_SUM = mpi4py_MPI.SUM if mpi4py_MPI is not None else None
    OPERATION_MIN = mpi4py_MPI.MIN if mpi4py_MPI is not None else None
    OPERATION_MAX = mpi4py_MPI.MAX if mpi4py_MPI is not None else None
    __instance = None

    def __new__(cls):
        if HighPerformanceComputingConfiguration.__instance is not None:
            return HighPerformanceComputingConfiguration.__instance

        HighPerformanceComputingConfiguration.__instance = super(HighPerformanceComputingConfiguration, cls).__new__(cls)
        return HighPerformanceComputingConfiguration.__instance

    def __init__(self):
        super(HighPerformanceComputingConfiguration, self).__init__()

        if mpi4py_mpi_comm_world is not None:
            self._mpi_num_procs = mpi4py_mpi_comm_world.Get_size()
        else:
            self._mpi_num_procs = 1

        if self._mpi_num_procs > 1:  # MPI parallel run
            self._is_parallel = True
            self._mpi_rank = mpi4py_mpi_comm_world.Get_rank()
        else:  # Sequential run
            self._is_parallel = False

    @property
    def is_parallel(self):
        return self._is_parallel

    @property
    def total_nb_procs(self):
        return self._mpi_num_procs

    @property
    def cpu_number(self):
        if self._is_parallel:
            return self._mpi_rank + 1  # ranges in [ 1 ; mpi_num_procs ]

        return 1

    def iter_cpu_number(self):
        for icpu in range(1, self._mpi_num_procs+1):
            yield icpu

    @property
    def process_root(self):
        return 0

    @property
    def is_root(self):
        if self._is_parallel:
            return self._mpi_rank == self.process_root
        return True

    @property
    def comm_world(self):
        if self._is_parallel:
            return mpi4py_mpi_comm_world
        return None


__all__ = ["HighPerformanceComputingConfiguration"]
