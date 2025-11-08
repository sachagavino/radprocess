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
:mod:`pymses.filters.cell_to_points` --- CellsToPoints filter module
--------------------------------------------------------------------

"""

import numpy
from pymses.core import Filter, DataSource, IsotropicExtPointDataset, PointDataset


class CellsToPoints(Filter):
    r"""
    AMR grid to cell list conversion filter
    Filters an AMR dataset and converts it into a point-based dataset

    source:
        AMR source
    include_nonactive_cells (default False):
        If True, the created PointDataset keeps non active cells (i.e. ghost cells)
    include_boundary_cells (default False):
        If True, boundary cells are included
    include_split_cells (default False):
        If True, the created PointDataset will include all points from
        intermediary AMR resolution level (i.e. cells that are refined).
        If False, only leaf cell values are converted (this save memory
        and computation time for cell_to_points splatting rendering)
    smallest_cell_level ``integer`` (default None):
        If not None, the cells that are too small (compared to this given
        level of resolution) are filtered.
    """

    def __init__(self, source, include_nonactive_cells=False, include_boundary_cells=False,
                 include_split_cells=False, smallest_cell_level=None):

        self.include_split_cells = include_split_cells
        self.include_boundary_cells = include_boundary_cells
        self.include_nonactive_cells = include_nonactive_cells
        self.smallest_cell_level = smallest_cell_level

        super(CellsToPoints, self).__init__(source)

    def source_type(self):
        return DataSource.PARTICLE_SOURCE

    def filtered_dset(self, dset):
        r"""
        Filters an AMR dataset and converts it into a point-based dataset

        Returns
        -------
        PointAbstractDataset source

        """
        ndim = dset.amr_header["ndim"]
        twotondim = 1 << ndim

        #############
        # Grid mask #
        #############
        if not self.include_nonactive_cells:
            # Skip inactive grids if requested
            grid_mask = dset.get_active_mask()
        else:
            if not self.include_boundary_cells:
                # Skip boundary grids if requested
                grid_mask = ~dset.get_boundary_mask()
            else:
                # Set the initial grid mask to True
                grid_mask = numpy.ones(dset.amr_struct["ngrids"], 'bool')
        if self.smallest_cell_level is not None:
            # Filter too small cells
            grid_mask *= dset.get_grid_levels() <= self.smallest_cell_level
        # Compute the grid levels, the cell centers and then the cell levels
        cell_centers = dset.get_cell_centers(grid_mask).reshape((-1, ndim))

        gl = dset.get_grid_levels(grid_mask)[:, numpy.newaxis]
        cell_sizes = numpy.repeat(0.5 ** gl, twotondim, axis=1).reshape((-1,))

        # Initialize the output PointDataset + cell size field
        pts = IsotropicExtPointDataset(cell_centers, cell_sizes)

        # Transfer the data arrays into the points
        for name in dset.scalars:
            data = dset[name]
            pts.add_scalars(name, data[grid_mask, :].reshape((-1,)))
        for name in dset.vectors:
            data = dset[name]
            pts.add_vectors(name, data[grid_mask, :].reshape((-1, data.shape[2])))
        for name in dset.multivalued:
            data = dset[name]
            pts.add_multivalued(name, data[grid_mask, :].reshape((-1, data.shape[2])))

        # Keep only the leaf cells if required
        if not self.include_split_cells:
            if self.smallest_cell_level is None:
                smallest_cell_size = 0.5 ** dset.amr_struct["readlmax"]
            else:
                smallest_cell_size = 0.5 ** self.smallest_cell_level
            sons = dset.amr_struct["son_indices"][grid_mask, :].reshape((-1,))
            leaf = (sons < 0) + (cell_sizes == smallest_cell_size)
            pts = pts.filtered_by_mask(leaf)

        return pts


class SplitCells(Filter):
    r"""
    Create point-based data from cell-based data by splitting the cell-mass
    into uniformly-distributed particles

    """

    def __init__(self, source, info, particle_mass):
        self.d_unit = info["unit_density"]
        self.l_unit = info["unit_length"]
        self.m_unit = info["unit_mass"]
        cell_points = CellsToPoints(source)
        super(SplitCells, self).__init__(cell_points)
        self.part_mass = particle_mass

    def filtered_dset(self, dset):
        r"""
        Split cell filtering method

        Parameters
        ----------
        dset : AbstractDataset

        Returns
        -------
        fdset : AbstractDataset
            filtered dataset

        """
        dx = dset.get_sizes()
        ndim = dset.points.shape[1]
        volume = dx ** ndim
        factor = (self.d_unit * self.l_unit ** ndim).express(self.m_unit)
        mcell = dset["rho"] * volume * factor
        dsets = []
        for size in numpy.unique(dx):  # Level-by-level treatment of the AMR cells
            mask0 = (dx == size)
            mm = mcell[mask0]
            pts = dset.points[mask0]
            # Number of particles of (at least) 'self.part_mass' mass in each cell
            nsamp = numpy.ceil(mm / self.part_mass).astype('i')
            new_pts = numpy.repeat(pts, nsamp, axis=0)
            nsamps = numpy.repeat(nsamp, nsamp)
            # if 1 particle in the cell, keep it at the center.
            # Otherwise, sample uniformly-distributed positions in the cell
            mask = (nsamps > 1)
            nr = new_pts[mask].shape[0]
            new_pts[mask, :] = new_pts[mask, :] + size * numpy.random.uniform(low=-0.5, high=0.5, size=(nr, ndim))
            # individual particle mass value
            msamp = mm / nsamp
            msamps = numpy.repeat(msamp, nsamp)
            # Create new dataset
            ds = PointDataset(new_pts)
            ds.add_scalars("mass", msamps)
            # individual particle velocity
            if "vel" in dset._all_data:
                vel = dset["vel"]
                vels = numpy.repeat(vel[mask0], nsamp, axis=0)
                ds.add_vectors("vel", vels)
            # individual particle Pressure energy
            if "P" in dset._all_data:
                Pr = dset["P"]
                vol = volume[mask0]
                volsamp = vol / nsamp
                volsamps = numpy.repeat(volsamp, nsamp)
                Prs = numpy.repeat(Pr[mask0] * volsamps, nsamp)
                ds.add_scalars("P", Prs)
            dsets.append(ds)
        new_dset = dsets[0].concatenate(dsets)
        return new_dset


__all__ = ["CellsToPoints", "SplitCells"]
