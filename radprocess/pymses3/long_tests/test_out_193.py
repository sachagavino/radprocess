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

import numpy
import os
import h5py

from pymses import RamsesOutput
from pymses.analysis.cube3d import CubeExtractor
from pymses.core import PointDataset
from pymses.filters import RegionFilter, PointFunctionFilter
from pymses.analysis import bin_cylindrical, bin_spherical
from pymses.analysis.point_sampling import PointSamplingProcessor
from pymses.analysis import Camera, DataMap, ScalarOperator, MaxLevelOperator, MinimumTemperatureOperator, \
    MassWeightedDensityOperator
from pymses.analysis import splatting, raytracing
from pymses.analysis.slicing import SliceMap
from pymses.utils.regions import Cylinder, Sphere
from pymses.utils import constants as C


class TestAquarius193(object):
    _datadir = None
    _pwd = None
    _iout = None
    _ro = None

    @classmethod
    def setup_class(cls):
        cls._datadir = os.path.join(os.path.expanduser("~"), "data/Aquarius/output")
        cls._iout = 193
        cls._ro = RamsesOutput(cls._datadir, cls._iout)
        cls._ro.define_amr_scalar_field("hydro", "rho", 0)
        cls._ro.define_amr_scalar_field("hydro", "P", 4)
        cls._pwd = os.path.abspath(os.path.dirname(__file__))

    def test_sample_points(self):
        # Galactic cylinder parameters
        gal_center = [0.567811, 0.586055, 0.559156]  # in box units
        gal_radius = 0.00024132905460547268  # in box units
        gal_thickn = 0.00010238202316595811  # in box units
        gal_normal = [-0.172935, 0.977948, -0.117099]  # Norm = 1


        # Prepare to read the density field only
        source = self._ro.amr_source(["rho"])

        # Cylinder region
        cyl = Cylinder(gal_center, gal_normal, gal_radius, gal_thickn)

        # AMR density field point sampling
        psp = PointSamplingProcessor(source)
        # psp.disable_multiprocessing()
        h5fname = os.path.join(self._pwd, "sample_points_cyl.h5")
        pdset = PointDataset.from_HDF5(h5fname)

        point_dset = psp.process(pdset.points, add_level=True, add_cell_center=True, interpolation=True)

        assert numpy.allclose(point_dset["rho"], pdset["rho"], rtol=1.0e-6)
        assert numpy.allclose(point_dset["cell_center"], pdset["cell_center"], rtol=1.0e-6)
        assert numpy.allclose(point_dset["level"], pdset["level"], rtol=1.0e-6)

    def test_cyl_profile(self):
        # Galactic cylinder parameters
        gal_center = [0.567811, 0.586055, 0.559156]  # in box units
        gal_radius = 0.00024132905460547268  # in box units
        gal_thickn = 0.00010238202316595811  # in box units
        gal_normal = [-0.172935, 0.977948, -0.117099]  # Norm = 1

        # Prepare to read the density field only
        source = self._ro.amr_source(["rho"])

        # Cylinder region
        cyl = Cylinder(gal_center, gal_normal, gal_radius, gal_thickn)

        # AMR density field point sampling
        numpy.random.seed(1652336)
        points = cyl.random_points(1.0e6)  # 1M sampling points
        psp = PointSamplingProcessor(source)
        # psp.disable_multiprocessing()
        point_dset = psp.process(points)

        def rho_weight_func(dset):
            return dset["rho"]

        r_bins = numpy.linspace(0.0, gal_radius, 200)

        # Profile computation
        rho_profile = bin_cylindrical(point_dset, gal_center, gal_normal, rho_weight_func, r_bins, divide_by_counts=True)

        # Plot
        # Geometrical midpoint of the bins
        length = self._ro.info["unit_length"].express(C.kpc)
        bins_centers = length * (r_bins[1:] + r_bins[:-1]) / 2.
        dens = self._ro.info["unit_density"].express(C.H_cc)

        # h5f = h5py.File("./long_tests/cyl_profile.h5", 'w')
        # h5f.create_dataset("cyl_profile", data=rho_profile)
        # h5f.close()

        h5f = h5py.File(os.path.join(self._pwd, "cyl_profile.h5"), 'r')
        rho_profileB = h5f["/cyl_profile"][...]
        h5f.close()

        # print rho_profile
        assert numpy.sum(rho_profile - rho_profileB) < 10e-6

    def test_sph_profile(self):
        # Halo parameters
        halo_center = [0.567811, 0.586055, 0.559156]  # in box units
        halo_radius = 0.00075  # in box units

        # Prepare to read the mass/epoch fields only
        source = self._ro.particle_source(["mass", "epoch"])

        # Sphere region
        sph = Sphere(halo_center, halo_radius)

        # Filtering particles
        point_dset = RegionFilter(sph, source)
        dm_filter = lambda dset: dset["epoch"] == 0.0
        dm_parts = PointFunctionFilter(dm_filter, point_dset)

        # Profile computation
        m_weight_func = lambda dset: dset["mass"]
        r_bins = numpy.linspace(0.0, halo_radius, 200)

        # Mass profile
        # This triggers the actual reading of the particle data files from disk.
        mass_profile = bin_spherical(dm_parts, halo_center, m_weight_func, r_bins, divide_by_counts=False)

        # Density profile
        sph_vol = 4.0 / 3.0 * numpy.pi * r_bins ** 3
        shell_vol = numpy.diff(sph_vol)
        rho_profile = mass_profile / shell_vol

        # Plot
        # Geometrical midpoint of the bins
        length = self._ro.info["unit_length"].express(C.kpc)
        bins_centers = (r_bins[1:] + r_bins[:-1]) / 2. * length
        dens = self._ro.info["unit_density"].express(C.Msun / C.kpc ** 3)

        # h5f = h5py.File("./long_tests/sph_profile.h5", 'w')
        # h5f.create_dataset("sph_profile", data=rho_profile)
        # h5f.close()

        h5f = h5py.File(os.path.join(self._pwd, "sph_profile.h5"), 'r')
        rho_profileB = h5f["/sph_profile"][...]
        h5f.close()

        # print rho_profile
        assert numpy.allclose(rho_profile, rho_profileB, rtol=1.0e-6)

    def test_slice_density(self):
        # AMR data source
        amr = self._ro.amr_source(["rho"])

        # Defining a Camera object
        cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', up_vector='y', region_size=[1., 1.],
                     map_max_size=100, log_sensitive=True)

        # Density field access operator
        rho_op = ScalarOperator(lambda dset: dset["rho"], self._ro.info["unit_density"])

        # Slice map computation
        datamap = SliceMap(amr, cam, rho_op, z=0.4)
        # create a density slice map at z=0.4 depth position

        # h5f = h5py.File("./long_tests/slice_density.h5", mode='w')
        # h5f.create_dataset("map", data=map)
        # h5f.close()

        h5f = h5py.File(os.path.join(self._pwd, "slice_density.h5"), 'r')
        mapB = h5f["/map"][...]
        h5f.close()

        # print map
        assert numpy.allclose(numpy.log10(datamap.map), mapB, rtol=1.0e-6)

    def test_rt_Tmin(self):
        # Map operator : minimum temperature along line-of-sight
        scal_op = MinimumTemperatureOperator("P", "rho", self._ro.info["unit_temperature"])

        # Map region
        center = [0.567111, 0.586555, 0.559156]
        axes = {"los": "z"}

        axname = "los"
        axis = "z"
        cam = Camera(center=center, line_of_sight_axis=axis, up_vector="y", region_size=[3.0E-3, 3.0E-3],
                     distance=1.5E-3, far_cut_depth=1.5E-3, map_max_size=100)

        # Map processing
        rt = raytracing.RayTracer(self._ro.amr_source(["rho", "P"]), self._ro.info, scal_op)
        dmap = rt.process(cam)


        Tmin_dmap_fname = os.path.join(self._pwd, "rt_min.h5")
        # dmap.save_HDF5(Tmin_dmap_fname)
        test_dmap = DataMap.from_HDF5(Tmin_dmap_fname)

        # from matplotlib import pyplot as P
        # P.figure()
        # P.imshow(numpy.log10(dmap.map.T), origin='lower')
        # P.colorbar()
        # P.figure()
        # P.imshow(numpy.log10(test_dmap.T), origin='lower')
        # P.colorbar()
        # P.show()
        assert dmap == test_dmap

    def test_rt_rho(self):
        # Map operator : mass-weighted density map
        scal_op = MassWeightedDensityOperator("rho", self._ro.info["unit_density"])

        # Map region
        center = [0.567811, 0.586055, 0.559156]
        axname = "los"
        axis = [-0.172935, 0.977948, -0.117099]
        cam = Camera(center=center, line_of_sight_axis=axis, up_vector="z", region_size=[3.0E-2, 3.0E-2],
                     distance=2.0E-2, far_cut_depth=2.0E-2, map_max_size=100)

        # Map processing
        rt = raytracing.RayTracer(self._ro.amr_source(["rho"]), self._ro.info, scal_op)
        dmap = rt.process(cam, surf_qty=False)

        # h5f = h5py.File("./long_tests/rt_rho.h5", 'w')
        # h5f.create_dataset("map", data=map)
        # h5f.close()

        h5f = h5py.File(os.path.join(self._pwd, "rt_rho.h5"), 'r')
        mapB = h5f["/map"][...]
        h5f.close()

        # print map
        assert numpy.allclose(dmap.map, mapB, rtol=1.0e-6)

    def test_rt_optical_depth(self):
        # Map region
        center = [0.567811, 0.586055, 0.559156]
        axname = "los"
        axis = [-0.172935, 0.977948, -0.117099]
        cam = Camera(center=center, line_of_sight_axis=axis, up_vector="z", region_size=[3.0E-2, 3.0E-2],
                     distance=2.0E-2, far_cut_depth=2.0E-2, map_max_size=100)

        # Map processing
        rt = raytracing.OpticalDepthTracer(self._ro.amr_source(["rho"]), self._ro.info, "rho")
        cells = rt.process(cam)

        tau = cells["tau"]
        from pymses.analysis.splatting.map_bin2d import histo2D
        pts = numpy.zeros_like(cells.points)
        pts[:, :2] , pts[:, 2] = cam.project_points(cells.points)
        centered_map_box = cam.get_map_box()
        map_range = numpy.array([[centered_map_box.min_coords[0], centered_map_box.max_coords[0]],
                                 [centered_map_box.min_coords[1], centered_map_box.max_coords[1]],
                                 [centered_map_box.min_coords[2], centered_map_box.max_coords[2]]])

        map = histo2D(pts, [256, 256], map_range, {"tau": tau})

        map_unit = self._ro.info['unit_density']*self._ro.info['unit_length']
        map = numpy.log10(map["tau"]*map_unit.express(C.Msun / C.kpc**2))

        # from matplotlib import pyplot as P
        # vmin = numpy.min(map[map > 0.])
        # P.imshow(map, origin='lower')
        # fo = FormatStrFormatter("$10^{%d}$")
        # offset = numpy.ceil(vmin) - vmin
        # lo = IndexLocator(1.0, offset)
        # cb = P.colorbar(ticks=lo, format=fo)
        # # Set colorbar lable
        # cb.set_label("$M_{\odot}/kpc^{2}$")
        # P.show()

        map_center = numpy.array([[6.74049499, 6.68199057, 7.06004931, 6.99901891, 6.7961903, 6.56412166],
                                  [6.71816739, 6.84314875, 7.04456357, 7.33523869, 6.68610451, 6.75278867],
                                  [6.75950766, 6.99763695, 7.36832528, 7.42688259, 7.20266501, 6.66591754],
                                  [7.00923344, 7.1890434,  7.41326778, 7.94251749, 7.41687803, 6.69849191],
                                  [6.86131479, 6.82153632, 7.34913798, 7.56796314, 7.43505078, 6.59677941],
                                  [6.46090957, 6.61134082, 6.6465145,  6.83919183, 6.53015767, 6.74407785]])

        assert numpy.allclose(map[125:131, 125:131], map_center, rtol=1.0e-6)

    def test_rt_lmax(self):
        # Map operator : max. AMR level of refinement along the line-of-sight
        scal_op = MaxLevelOperator()

        # Map region
        center = [0.567811, 0.586055, 0.559156]
        axname = "los"
        axis = [-0.172935, 0.977948, -0.117099]
        cam = Camera(center=center, line_of_sight_axis=axis, up_vector="z", region_size=[3.0E-2, 3.0E-2],
                     distance=2.0E-2, far_cut_depth=2.0E-2, map_max_size=100)

        # Map processing
        rt = raytracing.RayTracer(self._ro.amr_source(["rho"]), self._ro.info, scal_op)
        dmap = rt.process(cam)

        lvlmax_dmap_fname = os.path.join(self._pwd, "rt_lmax.h5")
        # dmap.save_HDF5(lvlmax_dmap_fname)
        test_lvlmax_dmap = DataMap.from_HDF5(lvlmax_dmap_fname)

        # from matplotlib import pyplot as P
        # P.figure()
        # P.imshow(dmap.map.T, origin='lower')
        # P.colorbar()
        # P.figure()
        # P.imshow(test_lvlmax_dmap.map.T, origin='lower')
        # P.colorbar()
        # P.show()
        assert dmap == test_lvlmax_dmap

    def test_fft_part(self):
        parts = self._ro.particle_source(["mass", "level"])

        # Map operator : mass
        scal_func = ScalarOperator(lambda dset: dset["mass"], self._ro.info["unit_mass"])

        # Map region
        center = [0.567811, 0.586055, 0.559156]

        axname = "los"
        axis = [-0.172935, 0.977948, -0.117099]

        cam = Camera(center=center, line_of_sight_axis=axis, up_vector="z", region_size=[5.0E-1, 4.5E-1],
                     distance=2.0E-1, far_cut_depth=2.0E-1, map_max_size=100)

        # Map processing
        mp = splatting.SplatterProcessor(parts, self._ro.info, scal_func)
        map = mp.process(cam, surf_qty=True)

        # h5f = h5py.File("./long_tests/fft_part.h5", 'w')
        # h5f.create_dataset("map", data=map)
        # h5f.close()

        # print map.map[25:30, 50:55]
        h5f = h5py.File(os.path.join(self._pwd, "fft_part.h5"), 'r')
        mapB = h5f["/map"][...]
        # print mapB[25:30, 50:55]
        # from matplotlib import pyplot as P
        # P.imshow(numpy.log10(map.data), cmap='Blues')
        # P.figure()
        # P.imshow(numpy.log10(mapB), cmap='Blues')
        # P.show()
        h5f.close()

        assert (map.map - mapB).all() < 10e-6

    def test_fft_amr(self):
        amr = self._ro.amr_source(["rho"])

        # Map operator : mass-weighted density map
        scal_func = MassWeightedDensityOperator("rho", self._ro.info["unit_density"])

        # Map region
        center = [0.567811, 0.586055, 0.559156]

        axname = "los"
        axis = [-0.172935, 0.977948, -0.117099]

        cam = Camera(center=center, line_of_sight_axis=axis, up_vector="z", region_size=[5.0E-1, 4.5E-1],
                     distance=2.0E-1, far_cut_depth=2.0E-1, map_max_size=100)

        # Map processing
        mp = splatting.SplatterProcessor(amr, self._ro.info, scal_func)
        dmap = mp.process(cam, surf_qty=True)

        # h5f = h5py.File("./long_tests/fft_amr.h5", 'w')
        # h5f.create_dataset("map", data=map)
        # h5f.close()

        h5f = h5py.File(os.path.join(self._pwd, "fft_amr.h5"), 'r')
        mapB = h5f["/map"][...]
        h5f.close()
        # print map.map[40:42, 50:52]
        # print mapB[40:42, 50:52]

        # from matplotlib import pyplot as P
        # P.imshow(numpy.log10(map.data), cmap='Blues')
        # P.figure()
        # P.imshow(numpy.log10(mapB), cmap='Blues')
        # P.show()

        assert numpy.allclose(dmap.map, mapB, rtol=1.0e-6)

    def test_amr2cube_cartesian(self):
        # backuped referenced data :
        cube_slice_disk = numpy.array([[3.31552056e+03, 1.81019991e+05, 1.44361632e+06, 1.06553898e+06, 6.67027973e+04],
                                       [2.77241677e+03, 2.23644459e+05, 1.96378876e+06, 1.13240799e+06, 2.02268839e+04],
                                       [2.65935222e+03, 2.91311908e+05, 2.69012497e+06, 9.79126647e+05, 1.11434823e+03],
                                       [3.44135190e+03, 4.00791164e+05, 3.06939634e+06, 7.41074921e+05, 2.28813715e+05],
                                       [4.71830648e+03, 5.20836967e+05, 2.41345555e+06, 5.34213344e+05, 5.71283161e+05]])

        gal_center = [0.567811, 0.586055, 0.559156]  # in box units
        gal_radius = 2.5e-4  # in box units
        gal_normal = [-0.172935, 0.977948, -0.117099]  # Norm = 1

        cam = Camera(center=gal_center, line_of_sight_axis=[0.5, 0.5, 0.5], up_vector="y", distance=2.5E-1,
                     size_unit=self._ro.info["unit_length"], region_size=[4.0E-1, 4.0E-1], far_cut_depth=2.5E-1,
                     map_max_size=256)

        source = self._ro.amr_source(["rho"])
        dens_op = ScalarOperator(lambda dset: dset["rho"], self._ro.info["unit_density"])
        ce = CubeExtractor(source, dens_op)
        dcube = ce.process(cam, cube_size=2.5*gal_radius, resolution=128)
        # dcube.save_fits("~/disk_cube.fits", axis_unit=C.kpc)
        # print dcube.data[60:65, 60:65, 64]
        assert numpy.allclose(dcube.data[60:65, 60:65, 64], cube_slice_disk, rtol=1.0e-6)

    @classmethod
    def teardown_class(cls):
        pass


__all__ = ["TestAquarius193"]
