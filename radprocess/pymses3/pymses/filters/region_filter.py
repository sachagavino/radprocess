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
:mod:`pymses.filters.region_filter` --- Region filter module
------------------------------------------------------------
"""
import pymses
from pymses.core import DataSource, SubsetFilter
from numpy import ceil


class RegionFilter(SubsetFilter):
    r"""
    Region Filter class. Filters the data contained in a given region
    of interest.

    Parameters
    ----------
    region : :class:`~pymses.utils.regions.Region`
        region of interest
    source : :class:`~pymses.core.sources.DataSource`
        data source

    """

    def __init__(self, region, source):
        self.region = region
        if source.dom_decomp is not None:
            map_list_full = source.dom_decomp.map_region(self.region)
            if pymses.hpc_config.is_parallel:  # Spread the workload on all available CPUs => split the domain list evenly
                # between MPI processes
                ndm = len(map_list_full)
                if ndm <= pymses.hpc_config.total_nb_procs:  # More MPI processes than domains
                    if pymses.hpc_config.cpu_number <= ndm:
                        map_list = [map_list_full[pymses.hpc_config.cpu_number-1]]
                    else:
                        map_list = []
                else:  # More domains than MPI processes
                    ndm_per_proc = int(ceil(ndm / float(pymses.hpc_config.total_nb_procs)))
                    ndomain_remain = ndm
                    npp = ndm_per_proc
                    j = 0
                    for iproc in pymses.hpc_config.iter_cpu_number():
                        if ndomain_remain / float(pymses.hpc_config.total_nb_procs - iproc + 1) == npp - 1:
                            npp -= 1
                        if iproc == pymses.hpc_config.cpu_number:  # MPI process domain subset
                            map_list = map_list_full[j:j + npp]
                            break
                        else:
                            j += npp
                            ndomain_remain -= npp
            else:  # Sequential run process all domains on the CPU
                map_list = map_list_full
        else:
            map_list = None
        super(RegionFilter, self).__init__(map_list, source)

    def filtered_dset(self, dset):
        if self.source.source_type() == DataSource.AMR_SOURCE:
            # update the active mask
            grid_mask = dset.get_active_mask()
            grid_center = dset.amr_struct["grid_centers"][grid_mask, :]
            filter_mask = self.region.contains(grid_center)
            dset.active_mask[grid_mask] = grid_mask[grid_mask] * filter_mask
            # TODO test corner points instead of grid centers
            return dset
        else:
            return dset.filtered_by_mask(self.region.contains(dset.points))


__all__ = ["RegionFilter"]
