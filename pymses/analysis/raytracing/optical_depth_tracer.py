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
r"""
:mod:`pymses.analysis.raytracing.optical_depth_tracer` --- Optical depth ray-tracing module
-------------------------------------------------------------------------------------------

"""
import numpy
from time import time

from pymses.core import IsotropicExtPointDataset
from pymses.filters import *
from ..data_processor import DataProcessor, DataProcess
from pymses.sources.ramses.octree import CameraOctreeDataset
from pymses.sources.ramses.domain_decomposition import HilbertDomainDecomp
from pymses.sources.ramses.tree_utils import tree_search
from ..operator import ScalarOperator
from ._raytrace import raytrace_amr



def optical_depth_raytracing_process_task(idomain, amr_source, ramses_info, required_level, operator, rays, verbose):
    # Get domain dataset
    dset = amr_source.get_domain_dset(idomain, verbose)

    # Fetch active mask & grid level arrays
    active_mask = dset.get_active_mask()
    g_levels = dset.get_grid_levels()

    # We do the processing only if needed, i.e. only if the amr level min of active cells in the dset is <= required_level
    if len(g_levels[active_mask]) > 0 and numpy.min(g_levels[active_mask]) <= required_level:
        ray_vectors, ray_origins, ray_length_max = rays
        cell_centers = dset.get_cell_centers()

        # Integer variables
        nrays= ray_origins.shape[0]
        ndim = dset.amr_header["ndim"]
        twotondim = (1 << ndim)
        nfunc = operator.nscal_func()
        ngrids = dset.amr_struct["ngrids"]
        max_read_level = required_level

        # Big AMR arrays
        igrids = numpy.empty(ngrids, 'i')
        icells = numpy.empty(ngrids, 'i1')
        iactive_mask = numpy.array(active_mask, 'i')
        sons = dset.amr_struct["son_indices"]

        # Integrand variables
        rho_data = numpy.empty((ngrids, twotondim), 'd')
        for (key, rho_func) in operator.iter_scalar_func():
            rho_data[...] = rho_func(dset)

        # Output arrays
        tau = numpy.zeros(nrays, 'd')
        ray_length_map = numpy.zeros(nrays, 'd')

        full_octree = True
        if ramses_info["dom_decomp"] is not None and isinstance(ramses_info["dom_decomp"], HilbertDomainDecomp):
            # iteration over the minimal grid description of the domain
            pos_blocks, order_list, noverlap = ramses_info["dom_decomp"].minimal_domain(dset.icpu)
            order = order_list.astype('i')
            nblocks = pos_blocks.shape[0] - noverlap
            if nblocks > 1 or order[0] != 0:
                search_dict = tree_search(dset.amr_struct, pos_blocks, order)
                igrids = search_dict["grid_indices"]
                icells = search_dict["cell_indices"]
                full_octree = False

        if full_octree:
            # Iteration on the full octree
            nblocks = 8
            igrids[0:8] = 0
            icells[0:8] = list(range(8))

        param_info_int = numpy.empty(7, 'i')
        param_info_int[0] = nrays
        param_info_int[1] = ndim
        param_info_int[2] = nblocks
        param_info_int[3] = max_read_level
        param_info_int[4] = 1
        param_info_int[5] = 0
        param_info_int[6] = 0

        raytrace_amr(tau, ray_length_map, ray_origins, ray_vectors, ray_length_max, cell_centers, sons, rho_data,
                     igrids, icells, iactive_mask, g_levels, param_info_int)

        return tau, ray_length_map

    return None, None


class OpticalDepthTracingProcess(DataProcess):
    def __init__(self, queues, source, ramses_info, operator, rays, req_level, verbose):
        tasks_queue, results_queue = queues
        super(OpticalDepthTracingProcess, self).__init__(tasks_queue, results_queue)
        self.source = source
        self.ramses_info = ramses_info
        self.operator = operator
        self.map_dict = {}
        self.required_level = req_level

        # Get rays info
        self.rays = rays
        ncells = self.rays[0].shape[0]
        self.tau = numpy.zeros(ncells, dtype='d')
        self.ray_length_cells = numpy.zeros(ncells, dtype='d')
        self.verbose = verbose

    def process_task(self, idomain):
        task_res = optical_depth_raytracing_process_task(idomain, self.source, self.ramses_info, self.required_level,
                                                         self.operator, self.rays, self.verbose)

        if task_res[0] is not None:
            self.tau += task_res[0]
        if task_res[1] is not None:
            self.ray_length_cells += task_res[1]

        return None

    @property
    def result(self):
        res = (self.tau, self.ray_length_cells)
        return res


class OpticalDepthTracer(DataProcessor):
    r"""
    Optical depth tracer processing class

    Parameters
    ----------
    source: : class:`~pymses.sources.ramses.sources.RamsesAmrSource`
        AMR data source.
    ramses_output_info: ``dict``
        RamsesOutput info dict.
    density_field: ``string``
        Density field name.
    verbose: ``bool``
        verbosity boolean flag. Default None.
    """
    def __init__(self, source, ramses_output_info, density_field, verbose=None):
        if not isinstance(density_field, str):
            raise AttributeError("'density_field' must be a valid field name (string).")
        op = ScalarOperator(lambda dset: dset[density_field], ramses_output_info["unit_density"])
        super(OpticalDepthTracer, self).__init__(source, op, amr_mandatory=True, verbose=verbose)
        self._ro_info = ramses_output_info
        self._cells_source = None
        self._required_level = None
        self._camera = None
        self._rays = None

    def process_factory(self, tasks_queue, results_queue, verbosity):
        return OpticalDepthTracingProcess((tasks_queue, results_queue), self._filtered_source, self._ro_info,
                                          self._operator, self._rays, self._required_level, verbosity)

    def process(self, camera, verbose=True):
        r"""
        Map processing method : ray-tracing through data cube

        Parameters
        ----------
        camera          : :class:`~pymses.analysis.camera.Camera`
            camera containing all the view params
        verbose		: ``boolean`` (default False)
            show more console printouts

        Returns
        -------
        TODO
        """
        begin_time = time()

        self._camera = camera

        # AMR DataSource preparation
        self._source.set_read_levelmax(self._camera.get_required_resolution())
        self._required_level = self._source.get_read_levelmax()

        # Data bounding box
        domain_bounding_box = camera.get_bounding_box()

        # Extended domain bounding box for big octree cells
        ext = 0.5 ** (self._ro_info["levelmin"])
        domain_bounding_box.min_coords = numpy.amax([domain_bounding_box.min_coords - ext, [0., 0., 0.]], axis=0)
        domain_bounding_box.max_coords = numpy.amin([domain_bounding_box.max_coords + ext, [1., 1., 1.]], axis=0)

        # Data spatial filtering
        self._filtered_source = RegionFilter(domain_bounding_box, self._source)
        cells_source = CellsToPoints(self._filtered_source).flatten(verbose)
        ncells = cells_source.npoints
        # Get rays info
        self._rays = self._camera.get_rays(custom_origins=cells_source.points)
        ntasks = self._filtered_source.ndataset

        # Maps initialisation
        tau = numpy.zeros(ncells, dtype='d')
        ray_length_cells = numpy.zeros(ncells, dtype='d')

        if isinstance(self._source, CameraOctreeDataset):
            # In this case there is only one dataset to process so:
            self._use_hilbert_domain_decomp = False
            self.disable_multiprocessing()

        if self.use_multiprocessing and ntasks >= 2:
            self.init_multiprocessing_run(self._filtered_source.ndataset)

            # Loop over domains : enqueue jobs
            for idomain in self._filtered_source.iter_idata():
                self.enqueue_task(idomain)

            for res in self.get_results():
                if res[0] is not None:
                    tau += res[0]
                if res[1] is not None:
                    ray_length_cells += res[1]
        else:
            for idomain in self._filtered_source.iter_idata():
                res = optical_depth_raytracing_process_task(idomain, self._filtered_source, self._ro_info,
                                                            self._required_level, self._operator, self._rays,
                                                            self._verbose)

                if res[0] is not None:
                    tau += res[0]
                if res[1] is not None:
                    ray_length_cells += res[1]

        print("Optical depth computation time = %.3fs" % (time() - begin_time))

        ray_length = self._rays[2]
        if not numpy.allclose(ray_length_cells, ray_length, rtol=1.0e-4):
            print("Warning : calculated ray lengths during the ray-tracing process are not always equal." \
                  "Some cells may be missing.")

        # Build point dataset with tau and ray length scalar fields
        pts = IsotropicExtPointDataset(cells_source.points, cells_source.get_sizes())
        pts.add_scalars("tau", tau)
        pts.add_scalars("length", ray_length_cells)

        return pts


__all__ = ["OpticalDepthTracer"]
