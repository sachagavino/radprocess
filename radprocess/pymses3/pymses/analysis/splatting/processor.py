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

import numpy as N
from time import time

from ..data_processor import DataProcessor
try:
    from ..data_processor import DataProcess
except ImportError:
    DataProcess = object

from pymses.filters import *
from pymses.core import DataSource, IsotropicExtPointDataset, PointDataset
from .convolution_kernels import GaussSplatter2DKernel
from .camera import ExtendedCamera, CameraFilter
from .map_bin2d import histo2D
from pymses.sources.ramses.octree import CameraOctreeDataset
from pymses.analysis import DataMap
from ..operator import MultiFieldOperator
from pymses.utils.constants import none as None_unit


def splatter_binning_process_task(idomain, source, convol_kernel, camera, operator, project_point_kwargs, verbose):
    # Get domain dataset
    dset = source.get_domain_dset(idomain, verbose)

    if dset.npoints == 0:  # No point in dset => fetch next dataset
        if verbose is not None and verbose:
            print("No interesting point in this dataset")
        return None, None, None
    #print "npoints = %i" % dset.npoints

    # Sizes of the points
    sizes = convol_kernel.get_size(dset)

    # Big cells => straightforward summation
    mask = (sizes > project_point_kwargs['big_cells_size'])
    if mask.any():
        big_cells_dset = dset.filtered_by_mask(mask)

        # Small cells => FFT processing
        mask = ~mask
        small_cells_dset = dset.filtered_by_mask(mask)
        if small_cells_dset.npoints == 0:
            return big_cells_dset, None, None
        sizes = convol_kernel.get_size(small_cells_dset)

        # Get projected (u,v,w) coordinates of the dataset points
        points = dset.points[mask]
        if project_point_kwargs["random_shift"]:
            points = PointDataset.add_random_shift(points)

        # => u/v/w Coordinates of the points
        pts, w = camera.project_points(points)
    else:
        big_cells_dset = None
        small_cells_dset = dset
        points = dset.points
        if dset.npoints == 0:
            return big_cells_dset, None, None

        if project_point_kwargs["random_shift"]:
            points = PointDataset.add_random_shift(points)

        # => u/v/w Coordinates of the points
        pts, w = camera.project_points(points)

    ################################ Weights fields ###########################
    # => change weights depending on depth, with a window function
    # centered_map_box = camera.get_map_box()
    # zmin = centered_map_box.min_coords[2]
    # zmax = centered_map_box.max_coords[2]
    # f = N.abs(w_coord - (zmin + zmax) / 2.) / (zmax - zmin) - 0.49
    # weights_fact = win_func(f, l=0.01)

    # => Operator-by-operator weight dict.
    weights = {}
    map_dict = {}
    cam_dict = {}
    for key, func in operator:
        weights[key] = func(small_cells_dset)# * weights_fact

    # Map 2D binning
    for size in N.unique(sizes):  # Level-by-level cell processing
        map_dict[size] = {}
        mask = (size == sizes)

        w = {}
        for key, func in operator:
            w[key] = weights[key][mask]

        # Get the camera
        r = N.ones((2, 3)) * size
        r[0, :] = -r[0, :]
        # Convolution kernel = gaussian type : extend the region by 3*sigma
        # to get 99% of the gaussian into the extended region (anti-periodic effect)
        if isinstance(convol_kernel, GaussSplatter2DKernel):
            r = 3. * r
        cam_dict[size] = ExtendedCamera(camera, r)
        c = cam_dict[size]

        # Camera view area
        centered_map_box = c.get_map_box()
        map_range = N.array([[centered_map_box.min_coords[0], centered_map_box.max_coords[0]],
                             [centered_map_box.min_coords[1], centered_map_box.max_coords[1]],
                             [centered_map_box.min_coords[2], centered_map_box.max_coords[2]]])

        # Camera map size
        nx_map, ny_map = c.get_map_size()

        # 2D Binning of the dataset points with a dict. of weights.
        map_dict[size] = histo2D(pts[mask, :], [nx_map, ny_map], map_range, w)

    res = (big_cells_dset, map_dict, cam_dict)
    return res


def _cumul_dict(source, target, fill_camera=False, map_sum=False):
    if map_sum:
        for var in source:
            map = source[var]
            if map is None:
                continue
            if var not in target:
                target[var] = map
            else:
                target[var] += map
        return

    for size in source:
        o = source[size]
        if size in target:
            if fill_camera:
                continue

            t = target[size]
            for var in o:
                map = o[var]
                if var not in t:
                    t[var] = map
                else:
                    t[var] += map
        else:
            target[size] = o
    return


class SplatterBinningProcess(DataProcess):
    def __init__(self, queues, source, convol_kernel, camera, operator, sample_point_kwargs, verbose):
        tasks_queue, results_queue = queues
        super(SplatterBinningProcess, self).__init__(tasks_queue, results_queue)
        self.source = source
        self.convol_kernel = convol_kernel
        self.camera = camera
        self.operator = operator
        self.cam_dict = {}
        self.sample_point_kwargs = sample_point_kwargs
        self.big_cells = []
        self.map_dict = {}
        self.verbose = verbose

    def process_task(self, idomain):
        task_res = splatter_binning_process_task(idomain, self.source, self.convol_kernel, self.camera, self.operator,
                                                 self.sample_point_kwargs, self.verbose)

        if task_res[0] is not None:
            self.big_cells.append(task_res[0])
        if task_res[1] is not None:
            _cumul_dict(task_res[1], self.map_dict)
        if task_res[2] is not None:
            _cumul_dict(task_res[2], self.cam_dict, fill_camera=True)

        return None

    @property
    def result(self):
        res = (self.big_cells, self.map_dict, self.cam_dict)
        return res


class SplatterFFTConvProcess(DataProcess):
    def __init__(self, queues, convol_kernel, camera_dict, verbose):
        tasks_queue, results_queue = queues
        super(SplatterFFTConvProcess, self).__init__(tasks_queue, results_queue)
        self.convol_kernel = convol_kernel
        self.camera_dict = camera_dict
        self.maps = {}
        self.verbose = verbose

    def process_task(self, item):
        size, map_dict = item
        cam = self.camera_dict[size]
        _cumul_dict(self.convol_kernel.convol_fft(size, map_dict, cam, verbose=self.verbose), self.maps, map_sum=True)
        return None

    @property
    def result(self):
        return self.maps


class SplatterProcessor(DataProcessor):
    r"""
    Splatting data processor class.

    Parameters
    ----------
    source: ``DataSource``
        data source
    ramses_output_info: ``dict``
        RamsesOutput info dict.
    op: :class:`~pymses.analysis.operator.AbstractOperator`
        physical scalar quantity data operator
    verbose: ``bool``
        verbosity boolean flag. Default None.
    """
    def __init__(self, source, ramses_output_info, op, verbose=None):
        super(SplatterProcessor, self).__init__(source, op, verbose=verbose)
        self._ro_info = ramses_output_info
        self._convol_kernel = None
        self._pre_flatten = False
        self._filtered_source = None
        self._camera = None
        self._project_point_kwargs = None
        self._use_camera_lvlmax = True
        self._cam_dict = None
        self._process_state = None

    def process_factory(self, tasks_queue, results_queue, verbosity):
        if self._process_state == 0:
            return SplatterBinningProcess((tasks_queue, results_queue), self._filtered_source, self._convol_kernel,
                                          self._camera, self._operator, self._project_point_kwargs, verbosity)
        elif self._process_state == 1:
            return SplatterFFTConvProcess((tasks_queue, results_queue), self._convol_kernel, self._cam_dict,
                                          verbosity)

    def _prepare_data(self, camera):
        """prepare data method : it computes the "self._filtered_source" source attribute
        for the process(...) method.
        The data are then filtered with the CameraFilter class

        Parameters
        ----------
        camera          : :class:`~pymses.analysis.Camera`
            camera containing all the view params, the filtering is done according to those param
        """
        ext_size = self._convol_kernel.get_max_size()
        # Make sure we have a points dset source
        if self._source.source_type() == DataSource.AMR_SOURCE:
            if isinstance(self._source, CameraOctreeDataset):
                rlev = min(camera.get_required_resolution(), self._ro_info["levelmax"])
                points_dset_source = CellsToPoints(self._source, smallest_cell_level=rlev)
            else:
                points_dset_source = CellsToPoints(self._source)
        elif not (isinstance(self._source, ExtendedPointFilter) or isinstance(self._source, IsotropicExtPointDataset)):
            points_dset_source = ExtendedPointFilter(self._source)
        else:
            points_dset_source = self._source

        # Max. AMR level to read
        cam_lvlmax = camera.get_required_resolution()
        if self._use_camera_lvlmax and cam_lvlmax < self._ro_info["levelmax"]:
            points_dset_source.set_read_levelmax(cam_lvlmax)
        else:
            points_dset_source.set_read_levelmax(self._ro_info["levelmax"])

        # Filter the points dset source according to the camera
        self._filtered_source = CameraFilter(points_dset_source, camera, ext_size,
                                             use_camera_lvlmax=self._use_camera_lvlmax)
        if self._pre_flatten:
            self._filtered_source = self._filtered_source.flatten(self._verbose)

    def process(self, camera, ker_conv=None, pre_flatten=False, use_camera_lvlmax=True, surf_qty=True,
                FFTkernelSizeFactor=1, random_shift=False):
        """Map processing method

        Parameters
        ----------
        camera          : :class:`~pymses.analysis.camera.Camera`
            camera containing all the view params
        ker_conv : :class:`~pymses.analysis.splatting.convolution_kernels.ConvolKernel'  (default None leads to use a
        GaussSplatterKernel) convolution kernel used for the map processing
        surf_qty        : ``boolean``
            whether the processed map is a surface physical quantity. If True, the map is divided by the surface of a
            camera pixel. Default True.
        FFTkernelSizeFactor  : ``int or float`` (default 1)
            allow to change the convolution kernel size by a multiply factor to adjust points size
        random_shift : ``boolean`` (default False)
            add a random shift to point positions to avoid seeing the grid on resulting image
        pre_flatten : ``boolean`` (default False)
            Option to flatten the data source (using multiprocessing if possible) before computing the map
            The filtered data are then saved into the "self._filtered_source" source attribute.
        use_camera_lvlmax : ``boolean`` (default True)
            Limit the transformation of the AMR grid to particles to AMR cells under the camera octree levelmax
            (so that visible cells are only the ones that have bigger size than the camera pixel size).

        Returns
        -------
        map : ``array``
            FFT-convolved processed map
        """
        # Default convolution kernel : 2D gaussian kernel
        if ker_conv is None:
            max_ker_size = 0.5 ** (self._ro_info["levelmin"] + 1)
            self._convol_kernel = GaussSplatter2DKernel(max_size=max_ker_size)
        else:
            self._convol_kernel = ker_conv

        self._pre_flatten = pre_flatten
        self._camera = camera
        self._filtered_source = None
        self._use_camera_lvlmax = use_camera_lvlmax
        self._process_state = 0

        self._convol_kernel.FFTkernelSizeFactor = FFTkernelSizeFactor

        self._prepare_data(camera)

        # Region size level
        region_level = camera.get_region_size_level()
        cs = 1. / 2 ** region_level

        # Map initial value
        nx, ny = camera.get_map_size()
        map = N.zeros((nx, ny))

        # Cell/point projection & 2D binning
        map_dict = {}
        maps = {}
        for key, func in self._operator:
            maps[key] = N.zeros((nx, ny))
        self._cam_dict = {}
        big_cells = []
        big_cells_size = cs / 4.

        self._project_point_kwargs = {"random_shift": random_shift}

        tInit0 = time()

        self._project_point_kwargs["big_cells_size"] = big_cells_size

        if self.use_multiprocessing and self._filtered_source.ndataset >= 2:
            # Multithreaded binning
            self.init_multiprocessing_run(self._filtered_source.ndataset)

            # Loop over domains : enqueue jobs
            for idomain in self._filtered_source.iter_idata():
                self.enqueue_task(idomain)

            for res in self.get_results():
                big_cells += res[0]
                _cumul_dict(res[1], map_dict)
                _cumul_dict(res[2], self._cam_dict, fill_camera=True)

        else:
            # Sequential binning
            for idomain in self._filtered_source.iter_idata():
                res = splatter_binning_process_task(idomain, self._filtered_source, self._convol_kernel, self._camera,
                                                    self._operator, self._project_point_kwargs, self._verbose)
                if res[0] is not None:
                    big_cells.append(res[0])
                if res[1] is not None:
                    _cumul_dict(res[1], map_dict)
                if res[1] is not None:
                    _cumul_dict(res[2], self._cam_dict, fill_camera=True)


        tInit1 = time()
        print("-> Level-by-level point/cell projection time = %.3fs" % (tInit1 - tInit0))

        # Now proceed to the FFT-convolution of the binned maps
        self._process_state += 1

        # FFT-convolution
        ntask = 0
        ker_sizes = []
        for size, md in map_dict.items():
            ker_sizes.append(size)
            ntask += len(md)
        ker_sizes.sort()
        ker_sizes.reverse()

        if self.use_multiprocessing and ntask >= 2:
            # Multiprocessing
            self.init_multiprocessing_run(ntask)

            # Loop over data map dict. : enqueue jobs
            for size in ker_sizes:
                self.enqueue_task((size, map_dict[size]))

            # Update maps dict with results from all processes
            for res in self.get_results():
                _cumul_dict(res, maps, map_sum=True)
        else:
            # Sequential FFT-convolution
            for size in ker_sizes:
                _cumul_dict(self._convol_kernel.convol_fft(size, map_dict[size], self._cam_dict[size],
                                                           verbose=self._verbose), maps,	map_sum=True)

        del map_dict
        del self._cam_dict

        tInit2 = time()
        print("-> FFT convolution of the binned maps : dt = %.3fs" % (tInit2 - tInit1))

        # Large points direct kernel summation
        if len(big_cells) > 0:
            dset_bcells = big_cells[0].concatenate(big_cells)
            print("-> Adding %d big cells/points coming from AMR cells bigger than : %g" %\
                  (dset_bcells.npoints, big_cells_size))

            # Map edges/center coordinates
            xc, yc = camera.get_pixels_coordinates_centers()

            # u/v/w Coordinates of the points
            pts, w = camera.project_points(dset_bcells.points)

            # Sizes of the points
            sizes = self._convol_kernel.get_size(dset_bcells)

            # change weights depending on depth, with a window function
            # centered_map_box = camera.get_map_box()
            # zmin = centered_map_box.min_coords[2]
            # zmax = centered_map_box.max_coords[2]
            # f = N.abs(w - (zmin + zmax) / 2.) / (zmax - zmin) - 0.49
            # weights_fact = win_func(f, l=0.01)

            for key, func in self._operator:
                weights = func(dset_bcells)
                for i in range(pts.shape[0]):
                    s = sizes[i]
                    m = self._convol_kernel.ker_func((xc - pts[i, 0]), (yc - pts[i, 1]), s)
                    maps[key] += m * weights[i] #* weights_fact[i]

            print("-> Big cells/points direct summation time = %.3fs" % (time() - tInit2))

        del self._filtered_source  # Free filtered data
        # Final Operator process
        nproc = True
        for m in maps.values():
            nproc *= (m == 0.0).all()
        if nproc:
            return map
        S = camera.get_pixel_surface()

        if isinstance(self._operator, MultiFieldOperator):
            # Build datamap dictionary with consistent units
            map_dict = self._operator.operation(maps)
            output_unit_dict = self._operator.output_unit

            datamap = None
            for key, map in map_dict.items():
                if key not in output_unit_dict:
                    print("Warning : undefined Unit for map key '%s'. Set to 'none' Unit." % key)
                    u = None_unit
                else:
                    u = output_unit_dict[key]
                if surf_qty:
                    m = map / S
                else:
                    m = map
                if datamap is None:
                    datamap = DataMap(m, camera, u, map_name=key)
                else:
                    datamap.add_scalar_map(m, key, u)
        else:  # Return single datamap with unit defined in the operator
            map = map + self._operator.operation(maps)
            del maps

            if surf_qty:
                map /= S
            datamap = DataMap(map, camera, self._operator.output_unit)

        return datamap


def win_func(x, l=0.1, c=9.0):
    """
    Windowing function which applies a continuous shift from 1.0 (x<= 0.0) to 0.0 (x>=l) for input values x.

    Parameters
    ----------
    x: ``numpy.ndarray``
        input coordinate array.
    l: ``float``
        size of the window. Default 0.1
    c: ``float``
        exponential curve stiffness. Default 9.0

    Returns
    -------
    y: ``numpy.ndarray``
        window function applied on the input `x`coordinate array.
    """
    y = N.zeros_like(x)
    y[(x <= 0.)] = 1.0

    # exp 1
    mask = ((x <= l / 2.) * (x > 0.))
    y[mask] = 0.5 + 0.5 * (1.0 - N.exp((x[mask] - l / 2.) / (l / (2. * c)))) / (1.0 - N.exp(-c))

    # exp 2
    mask = ((x <= l) * (x > l / 2.))
    y[mask] = 0.5 - 0.5 * (1.0 - N.exp(-(x[mask] - l / 2.) / (l / (2. * c)))) / (1.0 - N.exp(-c))

    return y


__all__ = ["SplatterProcessor"]
