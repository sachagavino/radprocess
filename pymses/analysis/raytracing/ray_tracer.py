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
:mod:`pymses.analysis.raytracing.ray_tracer` --- raytracing module
------------------------------------------------------------------

"""
import pymses
from ..camera import SphericalCamera
from ..datamap import DataMap
from ..operator import MultiFieldOperator
from pymses.filters import *
from pymses.utils import constants as C, Box, UnitCube
from ..data_processor import DataProcessor, DataProcess
from pymses.sources.ramses.octree import CameraOctreeDataset# , CameraOctreeDatasource
from pymses.sources.ramses.domain_decomposition import HilbertDomainDecomp
from pymses.sources.ramses.tree_utils import tree_search
from ._raytrace import raytrace_amr, raytrace_amr_spherical

import numpy
from time import time


def raytracing_process_task(idomain, amr_source, ramses_info, required_level, operator, rays, is_sph_proj, verbose):
    # Get domain dataset
    dset = amr_source.get_domain_dset(idomain, verbose)

    # Fetch active mask & grid level arrays
    active_mask = dset.get_active_mask()
    g_levels = dset.get_grid_levels()

    # We do the processing only if needed, i.e. only if the amr level min of active cells in the dset is <= required_level
    act_g_levels = g_levels[active_mask]
    if len(act_g_levels) > 0 and numpy.min(act_g_levels) <= required_level:
        if is_sph_proj:
            ray_vectors, ray_origins, ray_length_max, dS_sr, front_clip = rays
        else:
            ray_vectors, ray_origins, ray_length_max = rays
        cell_centers = dset.get_cell_centers()

        # Integer variables
        nrays= ray_origins.shape[0]
        ndim = dset.amr_header["ndim"]
        twotondim = (1 << ndim)
        nfunc = operator.nscal_func()
        ngrids = dset.amr_struct["ngrids"]

        # If operator is a 'maximum along the line-of-sight', process the AMR tree down the max. level of refinement
        if operator.is_max_alos():
            max_op = 1
            max_read_level = dset.amr_struct["readlmax"]
        else:
            max_op = 0
            if is_sph_proj:
                max_read_level = dset.amr_struct["readlmax"]
            else:
                max_read_level = required_level

        if operator.use_cell_dx():
            use_dx = 1
        else:
            use_dx = 0

        # Big AMR arrays
        igrids = numpy.empty(ngrids, 'i')
        icells = numpy.empty(ngrids, 'i1')
        iactive_mask = numpy.array(active_mask, 'i')
        sons = dset.amr_struct["son_indices"]

        # Integrand variables
        ifunc = 0
        scal_data = numpy.empty((ngrids, twotondim, nfunc), 'd')
        for (key, scal_func) in operator.iter_scalar_func():
            scal_data[:, :, ifunc] = scal_func(dset)
            ifunc += 1

        # Output arrays
        maps = numpy.zeros((nrays, nfunc), 'd')
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
        param_info_int[4] = nfunc
        param_info_int[5] = max_op
        param_info_int[6] = use_dx

        if not is_sph_proj:
            raytrace_amr(maps, ray_length_map, ray_origins, ray_vectors, ray_length_max, cell_centers, sons, scal_data,
                         igrids, icells, iactive_mask, g_levels, param_info_int)
        else:
            raytrace_amr_spherical(maps, ray_length_map, ray_origins, ray_vectors, ray_length_max, dS_sr, front_clip,
                                   cell_centers, sons, scal_data, igrids, icells, iactive_mask, g_levels, param_info_int)

        return maps, ray_length_map

    return None, None


class RayTracingProcess(DataProcess):
    def __init__(self, queues, source, ramses_info, camera, operator, rays, is_sph_proj, req_level, verbose):
        tasks_queue, results_queue = queues
        super(RayTracingProcess, self).__init__(tasks_queue, results_queue)
        self.source = source
        self.ramses_info = ramses_info
        self.camera = camera
        self.operator = operator
        self.max_op_alos = self.operator.is_max_alos()
        self.map_dict = {}
        self.is_sph_proj = is_sph_proj
        self.required_level = req_level

        # Get rays info
        self.rays = rays
        n_rays = self.rays[0].shape[0]
        # nx_map, ny_map = camera.get_map_size()
        nfunc = self.operator.nscal_func()
        self.maps = numpy.zeros((n_rays, nfunc), dtype='d')
        self.ray_length_map = numpy.zeros(n_rays, dtype='d')
        self.verbose = verbose

    def process_task(self, idomain):
        task_res = raytracing_process_task(idomain, self.source, self.ramses_info, self.required_level,
                                           self.operator, self.rays, self.is_sph_proj, self.verbose)

        if task_res[0] is not None:
            if self.max_op_alos:
                numpy.maximum(self.maps, task_res[0], self.maps)
            else:
                self.maps += task_res[0]
        if task_res[1] is not None:
            self.ray_length_map += task_res[1]

        return None

    @property
    def result(self):
        res = (self.maps, self.ray_length_map)
        return res


class RayTracer(DataProcessor):
    r"""
    RayTracer class

    Parameters
    ----------
    source: : class:`~pymses.sources.ramses.sources.RamsesAmrSource`
        AMR data source.
    ramses_output_info: ``dict``
        RamsesOutput info dict.
    op : :class:`~pymses.analysis.operator.AbstractOperator`
            physical quantity data operator.
    verbose: ``bool``
        verbosity boolean flag. Default None
    """
    def __init__(self, source, ramses_output_info, op, verbose=None):
        super(RayTracer, self).__init__(source, op, amr_mandatory=True, verbose=verbose)
        self._ro_info = ramses_output_info
        self._filtered_source = None
        self._camera = None
        self._rays = None
        self._is_spherical_projection = False
        self._required_level = None

    def process_factory(self, tasks_queue, results_queue, verbosity):
        return RayTracingProcess((tasks_queue, results_queue), self._filtered_source, self._ro_info,
                                 self._camera, self._operator, self._rays, self._is_spherical_projection,
                                 self._required_level, verbosity)

    def process(self, camera, surf_qty=True):
        r"""
        Map processing method : ray-tracing through Ramses octree datasets

        Parameters
        ----------
        camera : :class:`~pymses.analysis.camera.Camera`
            camera containing all the view params
        surf_qty : ``boolean``
            whether the processed map is a surface physical quantity. If False, the map is multiplied by the surface of
            a camera pixel. Default True.

        Returns
        -------
        dmap: :class:`~pymses.analysis.datamap.DataMap`
            ray-traced DataMap instance.
        """
        begin_time = time()

        self._camera = camera
        spherical = isinstance(camera, SphericalCamera)

        # AMR DataSource preparation
        if not self._operator.is_max_alos() and not spherical:
            self._source.set_read_levelmax(self._camera.get_required_resolution())
        self._required_level = self._source.get_read_levelmax()

        # Data bounding box
        domain_bounding_box = camera.get_bounding_box()

        # Get rays info
        if not spherical:
            self._rays = camera.get_rays()
            self._is_spherical_projection = False
        else:
            self._rays = camera.get_rays()
            self._is_spherical_projection = True

        n_rays = self._rays[1].shape[0]
        nx_map, ny_map = camera.get_map_size()

        # Extended domain bounding box for big octree cells
        ext = 0.5 ** (self._ro_info["levelmin"])
        domain_bounding_box.min_coords = numpy.amax([domain_bounding_box.min_coords - ext, [0., 0., 0.]], axis=0)
        domain_bounding_box.max_coords = numpy.amin([domain_bounding_box.max_coords + ext, [1., 1., 1.]], axis=0)

        # Data spatial filtering
        if domain_bounding_box == UnitCube and not pymses.hpc_config.is_parallel:
            self._filtered_source = self._source
        else:
            self._filtered_source = RegionFilter(domain_bounding_box, self._source)
        ntasks = self._filtered_source.ndataset

        # Maps initialisation
        maps = numpy.zeros((n_rays, self._operator.nscal_func()), dtype='d')
        ray_length_map = numpy.zeros(n_rays, dtype='d')

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
                    if self._operator.is_max_alos():
                        numpy.maximum(maps, res[0], maps)
                    else:
                        maps += res[0]
                if res[1] is not None:
                    ray_length_map += res[1]
        else:
            for idomain in self._filtered_source.iter_idata():
                res = raytracing_process_task(idomain, self._filtered_source, self._ro_info, self._required_level,
                                              self._operator, self._rays, self._is_spherical_projection, self._verbose)

                if res[0] is not None:
                    if self._operator.is_max_alos():
                        numpy.maximum(maps, res[0], maps)
                    else:
                        maps += res[0]
                if res[1] is not None:
                    ray_length_map += res[1]

        # ------------------------------------------ MPI process data reduction -------------------------------------- #
        if pymses.hpc_config.is_parallel:
            if pymses.hpc_config.is_root:
                reduced_maps = numpy.zeros_like(maps)
                reduced_ray_lengths = numpy.zeros_like(ray_length_map)
            else:
                reduced_maps = None
                reduced_ray_lengths = None
            if self._operator.is_max_alos():
                pymses.hpc_config.comm_world.Reduce(maps, reduced_maps, op=pymses.hpc_config.OPERATION_MAX,
                                                    root=pymses.hpc_config.process_root)
            else:
                pymses.hpc_config.comm_world.Reduce(maps, reduced_maps, op=pymses.hpc_config.OPERATION_SUM,
                                                    root=pymses.hpc_config.process_root)
            pymses.hpc_config.comm_world.Reduce(ray_length_map, reduced_ray_lengths, op=pymses.hpc_config.OPERATION_SUM,
                                                root=pymses.hpc_config.process_root)

            if pymses.hpc_config.is_root:
                maps[...] = reduced_maps[...]
                ray_length_map[...] = reduced_ray_lengths[...]
        # ------------------------------------------------------------------------------------------------------------ #

        # Only root CPU returns a DataMap (or unique CPU in sequential run, obviously)
        if not pymses.hpc_config.is_root:
            return None

        print("Ray trace process time = %.3fs" % (time() - begin_time))

        # Ray length validity check (consistent integrated ray lengths ?)
        ray_length = self._rays[2][:]
        if self._is_spherical_projection:
            msk = (self._rays[4] > 0.0) * (self._rays[4] <= ray_length)
            ray_length[msk] -= self._rays[4][msk]

        if not numpy.allclose(ray_length_map, ray_length, rtol=1.0e-4):
            print("!!! Warning !!! : calculated ray lengths during the ray-tracing process are not always equal." \
                  "Some cells may be missing.")

        ifunc = 0
        mapd = {}
        S = camera.get_pixel_surface()

        for key, func in self._operator.iter_scalar_func():
            if spherical:
                if self._operator.is_max_alos() or not surf_qty:
                    mapd[key] = maps[:, ifunc].reshape(nx_map, ny_map)
                else:
                    dS_sr = self._rays[3].reshape(nx_map, ny_map)
                    mapd[key] = maps[:, ifunc].reshape(nx_map, ny_map) / dS_sr
            else:
                if self._operator.is_max_alos() or surf_qty:
                    mapd[key] = maps[:, ifunc].reshape(nx_map, ny_map)
                else:
                    mapd[key] = S * maps[:, ifunc].reshape(nx_map, ny_map)
            ifunc += 1

        if isinstance(self._operator, MultiFieldOperator):
            # Build datamap dictionary with consistent units
            map_dict = self._operator.operation(mapd)
            output_unit_dict = self._operator.output_unit

            datamap = None
            for key, map in map_dict.iter_items():
                if key not in output_unit_dict:
                    print("Warning : undefined Unit for map key '%s'. Set to 'none' Unit." % key)
                    u = C.none
                else:
                    u = output_unit_dict[key]
                    if datamap is None:
                        datamap = DataMap(map, camera, u, map_name=key)
                    else:
                        datamap.add_scalar_map(map, key, u)
        else:  # Return single datamap with unit defined in the operator
            map = self._operator.operation(mapd)
            datamap = DataMap(map, camera, self._operator.output_unit)

        return datamap



# class OctreeRayTracer(DataProcessor):
# 	r"""
# 	RayTracerDir class
#
# 	Parameters
# 	----------
# 	ramses_output   : :class:`~pymses.sources.ramses.output.RamsesOutput`
# 		ramses output from which data will be read to compute the map
# 	field_list      : ``list`` of ``string``
# 		list of all the required AMR fields to read (see :meth:`~pymses.sources.ramses.output.RamsesOutput.amr_source`)
#
# 	"""
#
# 	def __init__(self, *args):
# 		nargs = len(args)
# 		assert nargs in [1, 2]
# 		if nargs == 1:
# 			dset = args[0]
# 			assert isinstance(dset, CameraOctreeDataset)
# 			super(OctreeRayTracer, self).__init__(dset)
# 		else:
# 			super(OctreeRayTracer, self).__init__(None)
# 			self.ro, self.field_list = args
#
# 	def process(self, op, camera, surf_qty=False, return_image=True, rgb=True, use_C_code=True):
# 		r"""
# 		Map processing method : directional ray-tracing through AMR tree
#
# 		Parameters
# 		op              : :class:`~pymses.analysis.operator.AbstractOperator`
# 			physical scalar quantity data operator
# 		camera          : :class:`~pymses.analysis.camera.Camera`
# 			camera containing all the view params
# 		surf_qty        : ``boolean`` (default False)
# 			whether the processed map is a surface physical quantity. If True, the map
# 			is divided by the surface of a camera pixel.
# 		return_image        : ``boolean`` (default True)
# 			if True, return a PIL image (when rgb option is also True), else it returns
# 			a numpy array map
# 		rgb        : ``boolean`` (default True)
# 			if True, this code use the camera.color_tf to compute a rgb image
# 			if False, this code doesn't use the camera.color_tf, and works like the
# 			standard RayTracer. Then it returns two maps : the requested map,
# 			and the AMR levelmax map
# 		use_C_code : ``boolean`` (default True)
# 			Our pure C code is faster than the (not well optimized) Cython code,
# 			and should give the same result
# 		"""
# 		if self.source is None:
# 			begin_time = time()
# 			# CameraOctreeDatasource creation
# 			source = self.ro.amr_source(self.field_list)
# 			esize = 0.5 ** (self.ro.info["levelmin"] + 1)
# 			cod = CameraOctreeDatasource(camera, esize, source, include_split_cells=False)
# 			self.source = cod.dset
# 			print "CameraOctreeDatasource loaded up to level", camera.get_required_resolution(), \
# 				"with ngrids =", self.source.amr_struct["ngrids"], \
# 				"(loading time = %.2fs" % (time() - begin_time), ")"
#
# 		# Get rays info
# 		t0 = time()
# 		nx_map, ny_map = camera.get_map_size()
# 		ray_vectors, ray_origins, ray_lengths = camera.get_rays()
# 		n_rays = ray_origins.shape[0]
# 		rlev = camera.get_required_resolution()
#
# 		begin_time = time()
# 		try:
# 			# try to use multiprocessing on rays
# 			##########################################################################
# 			# multiprocessing ray tracing, pixel distribution (sort first rendering) #
# 			##########################################################################
# 			from multiprocessing import Process, cpu_count, Pipe
# 			from pymses.utils import misc
#
# 			NUMBER_OF_PROCESSES = min(n_rays / 10000 + 1, cpu_count(), misc.NUMBER_OF_PROCESSES_LIMIT)
# 			if NUMBER_OF_PROCESSES == 1:
# 				raise (Exception)  # don't use multiprocessing
# 			s = n_rays / NUMBER_OF_PROCESSES + 1
#
# 			def process_dset(child_conn, i_1, i_2):
# 				# Utility method for multiprocessing ray tracing
# 				maps = numpy.zeros((s + 1, 3), dtype='d')
# 				maps = ray_trace_octree(self.source, ray_origins[i_1:i_2], ray_vectors[i_1:i_2], \
# 										ray_lengths[i_1:i_2], op, camera.color_tf, level_max=rlev, verbose=False,
# 										rgb=rgb, \
# 										use_C_code=use_C_code)
# 				child_conn.send(maps)
# 				child_conn.close()
#
# 			# Start worker processes
# 			parent_conn = []
# 			for i in range(NUMBER_OF_PROCESSES):
# 				p_c, child_c = Pipe()
# 				parent_conn.append(p_c)
# 				Process(target=process_dset, args=(child_c, i * s, (i + 1) * s)).start()
#
# 			# Get results
# 			I = numpy.zeros((n_rays, 3), dtype='d')
# 			for i in range(NUMBER_OF_PROCESSES):
# 				I[i * s:(i + 1) * s] = parent_conn[i].recv()
# 		except Exception:
# 			print 'No multiprocessing...'
# 			I = ray_trace_octree(self.source, ray_origins, ray_vectors, ray_lengths, \
# 								 op, camera.color_tf, rgb=rgb, level_max=rlev, use_C_code=use_C_code)
# 		print "Octree ray trace processing time = %.3fs" % (time() - begin_time)
# 		if rgb:
# 			shape = I.shape
# 			I = I.reshape(shape[0] * shape[1])
# 			if sum(I != I) != 0:
# 				print "Error : There are", sum(I != I), " NaN value"
# 			if not return_image:
# 				return I.reshape((nx_map, ny_map, 3))
# 			else:
# 				from PIL import Image as I
#
# 				map = I.reshape((nx_map * ny_map, 3))
# 				# map = (map - min(map)) / (max(map) - min(map))
# 				map[:, 0] = (map[:, 0] - min(map[:, 0])) / (max(map[:, 0]) - min(map[:, 0]))
# 				map[:, 1] = (map[:, 1] - min(map[:, 1])) / (max(map[:, 1]) - min(map[:, 1]))
# 				map[:, 2] = (map[:, 2] - min(map[:, 2])) / (max(map[:, 2]) - min(map[:, 2]))
# 				map = numpy.asarray(map * 255, dtype='i')
# 				R_band = I.new("L", (nx_map, ny_map))
# 				R_band.putdata(map[:, 0])
# 				# import pylab as P
# 				# P.imshow(map[:,2].reshape((nx_map,ny_map)))
# 				# P.show()
# 				G_band = I.new("L", (nx_map, ny_map))
# 				G_band.putdata(map[:, 1])
# 				B_band = I.new("L", (nx_map, ny_map))
# 				B_band.putdata(map[:, 2])
# 				img = I.merge("RGB", (R_band, G_band, B_band))
# 				return img.rotate(90)
# 		else:
# 			ifunc = 0
# 			map_dict = {}
# 			S = camera.get_pixel_surface()
# 			for key, func in op.iter_scalar_func():
# 				if (op.is_max_alos() + surf_qty):
# 					map_dict[key] = I[:, ifunc].reshape(nx_map, ny_map)
# 				else:
# 					map_dict[key] = S * I[:, ifunc].reshape(nx_map, ny_map)
# 				ifunc += 1
# 			map = op.operation(map_dict)
#
# 			levelmax_map = I[:, 2].reshape(nx_map, ny_map)
# 			return map, levelmax_map


__all__ = ["RayTracer"]  # , "OctreeRayTracer"]
