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
:mod:`pymses.analysis.ppv.processor` --- PPV datacube generator module
----------------------------------------------------------------------

"""
import numpy as N
from time import time

from pymses.sources.ramses.tree_utils import tree_search
from pymses.sources.ramses.domain_decomposition import HilbertDomainDecomp
from pymses.filters import RegionFilter
from .cube import PPVCube
from ..data_processor import DataProcessor
from ._ray_cast_ppv import ray_cast_ppv


# from utils import args, tools


class PPVProcessor(DataProcessor):
    """
    PPV datacube generator class

    Parameters
    ----------
    ramses_output: :class:`~pymses.sources.ramses.output.RamsesOutput`
        RAMSES simulation output on which the PPV datacube must be processed
    rho_fname: `string`
        AMR density field name
    vel_fname: `string`
        AMR velocity field name
    """

    def __init__(self, ramses_output, rho_fname, vel_fname):
        self._ro = ramses_output
        self._vel_fname = vel_fname
        self._rho_fname = rho_fname

        # AMR DataSource preparation
        super(PPVProcessor, self).__init__(self._ro.amr_source([self._rho_fname, self._vel_fname]), None)

    def process(self, camera, vrange=None, vunit=None, nv=100, sigma_unit=None, use_hilbert_domain_decomp=True):
        """
        PPV datacube processing method

        Paramters
        ---------
        camera: :class:`~pymses.analysis.camera.Camera`
            view camera
        vrange: `list` of `float`
            velocity range [vmin, vmax]. Default [-50.0, 50.0]
        vunit: :class:`~pymses.utils.constants.unit.Unit`
            velocity range unit. Default: km/s
        nv: `int`
            number of velocity bins (must be > 2). Default 100.
        sigma_unit: :class:`~pymses.utils.constants.unit.Unit`
            datacube values surface density unit. Default: mH/cm^2.
        use_hilbert_domain_decomp: `bool`
            Make use of the Hilbert domain decomposition to speed up the ray-casting algorithm ? Default: True.

        Returns
        -------
        cube: :class: `~pymses.analysis.ppv.cube.PPVCube`
            PPV datacube object
        """
        begin_time = time()

        # Datacube initialisation
        ppv_cube = PPVCube(camera, sigma_unit, vrange, nv, vunit)
        base_Sunit = self._ro.info["unit_density"] * self._ro.info["unit_length"]
        Scoeff = base_Sunit.express(ppv_cube.unit)
        vcoeff = ppv_cube.velocity_unit.express(self._ro.info["unit_velocity"])
        vmin = ppv_cube.min_velocity * vcoeff
        vmax = ppv_cube.max_velocity * vcoeff

        rlev = camera.get_required_resolution()
        self._source.set_read_levelmax(rlev)

        # Get rays info
        ray_vectors, ray_origins, ray_lengths = camera.get_rays()
        n_rays = ray_origins.shape[0]

        # Data spatial filtering
        domain_bounding_box = camera.get_bounding_box()
        # Extended domain bounding box for big octree cells
        ext = 0.5 ** (self._ro.info["levelmin"])
        domain_bounding_box.min_coords = N.amax([domain_bounding_box.min_coords - ext, [0., 0., 0.]], axis=0)
        domain_bounding_box.max_coords = N.amin([domain_bounding_box.max_coords + ext, [1., 1., 1.]], axis=0)
        rsource = RegionFilter(domain_bounding_box, self._source)

        ray_length_maps = N.zeros(n_rays, dtype='d', order='C')
        uaxis, vaxis, los_axis = camera.get_camera_axis()

        # TODO multiprocessing

        ##########################
        # Sequential ray tracing #
        ##########################
        for icpu in rsource.data_list:
            dset = rsource.get_domain_dset(icpu)
            active_mask = dset.get_active_mask()
            g_levels = dset.get_grid_levels()
            sons = dset.amr_struct["son_indices"]
            cell_centers = dset.get_cell_centers()

            use_default_blocks = True
            if self._ro.info["dom_decomp"] is not None:
                # Make good use of the Hilbert domain decomposition minimal cubic paving
                if use_hilbert_domain_decomp and isinstance(self._ro.info["dom_decomp"], HilbertDomainDecomp):
                    # iteration over the minimal grid description of the domain
                    pos_blocks, order_list, noverlap = self._ro.info["dom_decomp"].minimal_domain(dset.icpu)
                    nblocks = pos_blocks.shape[0] - noverlap
                    if nblocks > 1 or order_list[0] != 0:
                        search_dict = tree_search(dset.amr_struct, pos_blocks, order_list)
                        igrids = search_dict["grid_indices"]
                        icells = search_dict["cell_indices"]
                        use_default_blocks = False

            if use_default_blocks:
                # iteration on the full octree
                nblocks = 8
                igrids = N.zeros(nblocks, dtype='i')
                icells = N.arange(nblocks, dtype='i1')

            # We do the processing only if needed, i.e. only if the amr level min of active cells in the dset is <= rlev
            if len(g_levels[active_mask]) > 0 and N.min(g_levels[active_mask]) <= rlev:
                # Apply the sign correction to the velocity (positive velocities are for gas moving away from the observer)
                velocity = - N.sum(dset[self._vel_fname] * los_axis[N.newaxis, N.newaxis, :], axis=-1)
                # Apply the surface density unit conversion factor
                density = dset[self._rho_fname] * Scoeff
                # Cast active mask boolean numpy array into int array
                act_mask = active_mask.astype('i')
                ray_cast_ppv(ppv_cube.data, density, velocity, ray_length_maps, ray_origins, ray_vectors, ray_lengths,
                             vmin, vmax, cell_centers, sons, nblocks, igrids, icells, rlev, act_mask, g_levels)

        print("PPV computation time = %.3fs" % (time() - begin_time))

        different_ray_length = N.unique(ray_length_maps)
        equal = True
        if len(different_ray_length) > 1:
            if (N.max(different_ray_length) - N.min(different_ray_length)) > ray_lengths[0] * 10e-3:
                equal = False
                print("Calculated ray lengths during the ray trace process are not always equal")
                if len(different_ray_length) < 5:
                    print("ray_lengths[0] =", ray_lengths[0])
                    for value in different_ray_length:
                        print("There are", sum(ray_length_maps == value), "ray(s) with value", value)
        if equal:
            print("Calculated ray lengths during the ray trace process are all equal : visualized volume is complete.")

        return ppv_cube


__all__ = ["PPVProcessor"]
