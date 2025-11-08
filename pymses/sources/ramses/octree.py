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
:mod:`pymses.sources.ramses.octree` --- RAMSES octree module
------------------------------------------------------------
"""
import numpy
from time import time

from ._octree_utils import sample_octree_dataset
from . import tree_utils, filename_utils
from .amr import read_ramses_amr_file
from .amrdata import read_ramses_amr_data_file
from pymses.core.datasets import AbstractDataset, PointDataset
from pymses.core import DataSource
from pymses.core.reader import MultiFileDataReader
from pymses.filters import RegionFilter
from pymses.utils.point_utils import corner_points
from pymses.utils.regions import Box, Sphere


class OctreeDataset(AbstractDataset):
    """
    Octree dataset

    Parameters
    ----------
    amr_dicts: ``tuple`` (amr_header, amr_struct)
        tuple containing the amr_header and amr_struc dictionaries describing the octree structure
    """
    def __init__(self, amr_dicts):
        self.amr_header, self.amr_struct = amr_dicts
        super(OctreeDataset, self).__init__()
        self.active_mask = None
        self.icpu = 0

    def source_type(self):
        r"""
        Returns
        -------
        odset_type: DataSource.AMR_SOURCE

        """
        return DataSource.AMR_SOURCE

    def get_cell_centers(self, grid_mask=None):
        r"""
        Returns
        -------
        cell_centers : ``array``
            AMR cell center coordinates array

        """
        # Grid centers
        gc = self.amr_struct["grid_centers"]
        if grid_mask is not None:
            gc = gc[grid_mask, :]
        # Grid levels
        gl = self.get_grid_levels(grid_mask)

        # Compute the cell centers
        cell_centers = tree_utils.cell_centers_from_grid_centers(gc, gl)
        return cell_centers

    def get_grid_levels(self, grid_mask=None):
        r"""
        Returns
        -------

        g_levels : ``array``
            the grid levels array

        """
        if "cell_levels" in self.amr_struct:
            gl = self.amr_struct["cell_levels"]
        else:
            gl = tree_utils.grid_levels(self.amr_struct)
            self.amr_struct["cell_levels"] = gl
        if grid_mask is None:
            return gl
        else:
            return gl[grid_mask]

    def get_active_mask(self):
        r"""
        Returns
        -------
        mask = ``array`` of ``bool``
            Active grids mask

        """
        if self.active_mask is None:
            self.active_mask = numpy.ones(self.amr_struct["ngrids"], 'bool')
        return self.active_mask

    def sample_points(self, points, max_search_level=None, add_level=False, add_cell_center=False,
                      interpolation=False):
        r"""
        AMR grid point sampling method

        Parameters
        ----------
        points : ``array``
            sampling points coordinates array
        max_search_level : ``int`` or None (default None)
            max. AMR level to read during sampling
            (see :func:'~pymses.sources.ramses.tree_utils.tree_search' for details)
        add_level : ``boolean`` (default False)
            whether we need to add a `level` field in the returned dataset containing
            the value of the AMR level the sampling points fall into
        add_cell_center : ``boolean`` (default False)
            whether we need to add a `cell_center` field in the returned dataset containing
            the coordinates of the AMR cell center the sampling points fall into
        interpolation : ``boolean`` (default False)
            Experimental : attempt to compute bi/tri-linear interpolation

        Returns
        -------
        dset : ``PointDataset``
            point-based sampled values dataset of the available AMR fields

        """
        # Create an empty PointDataset based on 'points'
        pts = numpy.asarray(points)
        npoints = pts.shape[0]
        ndim = pts.shape[1]
        point_dset = PointDataset(pts)

        if max_search_level is None:
            search_levelmax = self.amr_struct["readlmax"]
        else:
            if not isinstance(max_search_level, int):
                raise AttributeError("'max_search_level'attribute is not a valid integer value.")
            if max_search_level > self.amr_struct["readlmax"]:
                search_levelmax = self.amr_struct["readlmax"]
            else:
                search_levelmax = max_search_level

        ngrids = self.amr_struct["ngrids"]
        fl = []
        for scal in self.scalars:
            fl.append(self[scal])
        for vect in self.vectors:
            for idim in range(self[vect].shape[2]):
                fl.append(self[vect][:, :, idim])
        for mv in self.multivalued:
            for ivar in range(self[mv].shape[2]):
                fl.append(self[mv][:, :, ivar])

        nscalars = len(fl)
        scalars_array = numpy.concatenate(fl)

        # Big AMR arrays
        cell_centers = self.get_cell_centers()
        grid_centers = self.amr_struct["grid_centers"]
        son_indices = self.amr_struct["son_indices"]

        extracted_data = sample_octree_dataset(pts, npoints, scalars_array, nscalars, grid_centers, cell_centers,
                                               son_indices, ngrids, ndim, search_levelmax, int(interpolation),
                                               int(add_level), int(add_cell_center))

        i = 0
        for scal in self.scalars:
            point_dset.add_scalars(scal, extracted_data[:, i])
            i += 1
        for vect in self.vectors:
            n2 = self[vect].shape[2]
            point_dset.add_vectors(vect, extracted_data[:, i:i + n2])
            i += n2
        for mv in self.multivalued:
            n2 = self[mv].shape[2]
            point_dset.add_multivalued(mv, extracted_data[:, i:i + n2])
            i += n2

        # Extract the value of the cell level in which the points are sampled
        if add_level:
            point_dset.add_scalars("level", extracted_data[:, i])
            i += 1

        # Extract the cel center coordinates in which the points are sampled
        if add_cell_center:
            point_dset.add_vectors("cell_center", extracted_data[:, i:i+ndim])

        return point_dset

    def _h5_serialize(self, h5group, **kwargs):
        """
        Serialize the octree dataset into a HDF5 file.

        Parameters
        ----------
        h5fg: ``h5py.Group``
            HDF5 group object where to store the octree dataset object.
        """
        # Save dataset icpu number
        h5group.attrs["icpu"] = self.icpu

        # Record amr_header and amr_struct dict contents
        for dict_name, d in [("amr_header", self.amr_header), ("amr_struct", self.amr_struct)]:
            group = h5group.create_group(dict_name)
            for data_name, amr_data in d.items():
                if isinstance(amr_data, numpy.ndarray):
                    group.create_dataset(data_name, data=amr_data)
                else:
                    group.attrs[data_name] = amr_data

        super(OctreeDataset, self)._h5_serialize(h5group, **kwargs)

    @classmethod
    def _h5_deserialize(cls, h5group, version):
        """
        Read a octree dataset from a HDF5 file.

        Parameters
        ----------
        h5group: ``h5py.Group``
            Main group to deserialize the OctreeDataset from.
        version: ``int``
            Version of the OctreeDataset class to deserialize.

        Returns
        -------
        dset: :class:`~pymses.core.datasets.OctreeDataset`
            new octree dataset instance
        """
        # Check PointDataset object version written in HDF5 file
        if version == cls._pymses_h5_version:  # Current PointDataset version (h5py syntax)
            amr_dicts = []
            for dict_name in ["amr_header", "amr_struct"]:
                d = {}
                group = h5group[dict_name]
                for data_name, amr_data in group.items():  # Read numpy data arrays
                    d[data_name] = amr_data[...]
                for attr_name in group.attrs:  # Read attributes
                    d[attr_name] = group.attrs[attr_name]

                amr_dicts.append(d)

            # Init Octree dataset
            odset = cls(tuple(amr_dicts))
            odset.set_read_levelmax(odset.amr_struct["readlmax"])

            # Set octree dataset icpu number
            odset.icpu = h5group.attrs['icpu']

            # Read fields
            odset._read_fields_from_HDF5(h5group)
        else:
            raise AttributeError("Unknown OctreeDataset version (%d)" % int(version))

        return odset


class RamsesOctreeDataset(OctreeDataset):
    r"""
    RAMSES octree dataset class

    contains all the relevant information about the AMR tree structure

    """
    def get_active_mask(self):
        r"""
        Returns
        -------
        mask = ``array`` of ``bool``
            Active grids mask

        """
        if self.active_mask is None:
            self.active_mask = numpy.zeros(self.amr_struct["ngrids"], 'bool')
            offset = 0
            icpu = self.icpu
            ngridlevel = self.amr_struct["ngridlevel"]
            ngridbound = self.amr_struct["ngridbound"]
            for ilevel in range(self.amr_struct["readlmax"]):
                nbefore = ngridlevel[:icpu, ilevel].sum()
                nkeep = ngridlevel[icpu, ilevel].sum()
                ngrids = ngridlevel[:, ilevel].sum() + ngridbound[:, ilevel].sum()
                self.active_mask[offset + nbefore:offset + nbefore + nkeep] = True
                offset += ngrids
        return self.active_mask

    def get_boundary_mask(self):
        r"""
        Returns
        -------
        mask = ``array`` of ``bool``
            Boundary grids mask

        """
        mask = numpy.zeros(self.amr_struct["ngrids"], 'bool')
        offset = 0
        ngridlevel = self.amr_struct["ngridlevel"]
        ngridbound = self.amr_struct["ngridbound"]
        for ilevel in range(self.amr_struct["readlmax"]):
            nskip = ngridlevel[:, ilevel].sum()
            ngrids = nskip + ngridbound[:, ilevel].sum()
            mask[offset + nskip:offset + ngrids] = True
            offset += ngrids
        return mask

    def get_idomain_grid(self):
        r"""
        Returns
        -------
        idomain_grid = ``array`` of ``int``
            idomain grids array : every amr oct in amr_struct match one idomain and cpu file.
            ! WARNING : icpu = idomain +1 So we have for icpu between 1 to nproc:
            dset = source.get_domain_dset(icpu)
            domain_grid = dset.get_idomain_grid()
            idomain = icpu - 1
            sum(domain_grid == idomain) == sum(dset.get_active_mask())

        """
        ngrids = self.amr_struct["ngrids"]
        ngridlevel = self.amr_struct["ngridlevel"]
        ngridbound = self.amr_struct["ngridbound"]
        idomain_grid = numpy.zeros(ngrids, "i")
        ngridarray = numpy.zeros((ngridlevel.shape[0] + ngridbound.shape[0], ngridlevel.shape[1]), "i")
        ngridarray[:ngridlevel.shape[0], :] = ngridlevel
        ngridarray[ngridlevel.shape[0]:, :] = ngridbound

        i = 0
        for cpu_lvl in ngridarray.transpose():
            for idomain, noct in enumerate(cpu_lvl):
                if noct != 0:
                    for xoct in range(i, i + noct):
                        idomain_grid[xoct] = idomain
                    i += noct

        assert i == ngrids
        return idomain_grid


class RamsesOctreeReader(MultiFileDataReader):
    r"""
    Ramses AMR octree reader class. Provides a single read_file() method to fetch the RamsesOctreeDataset corresponding
    to a Hilbert domain.

    Parameters
    ----------
    output_repos: ``string``
        Ramses simulation output directory path.
    iout: ``int``
        Ramses output number.
    ivars_descrs_by_file: ``dict``
        (ivars, field) file type-based dictionary.
    swap: ``bool``
        swap bytes if True (to match different endianness).
    grav_compat: ``bool``
        Old Ramses versions (prior to commit 5dd90f3, 2012-10-04) and new Ramses versions (posterior to commit bce4454,
        2015-07-09) work with nvar_file (number of scalar fields written in the file) parameter written in file header :
            \_,-> gravitational attraction vector field only (old version)
            \_,-> potential scalar field (phi) + gravitational attraction vector field (new version)

        For Ramses versions BETWEEN 2012-10-04 AND 2015-07-09, ndim + 1 variables were saved but only ndim was written
        in the file header.
        To read output files written by this flavor of Ramses, set the `grav_compat` attribute to True. (Default False).
    verbose: ``bool``
        verbosity flag. Default True.
    """

    def __init__(self, output_repos, iout, ivars_descrs_by_file, swap=False, grav_compat=False, verbose=True):

        self._output_repos = output_repos
        self._iout = iout
        self.ivars_descrs_by_file = ivars_descrs_by_file
        self._swap = swap
        self._grav_compat = grav_compat
        self._verbose = verbose

    def read_file(self, icpu, read_lmax, verbose=None):
        """Read the particles data from the .part file and return a PointDataset
        containing the  user-defined fields
        Parameters
        ----------
        icpu; ``int``
            cpu index of the particle files to read
        read_lmax              : :  ``boolean`` (default None)
            level max to read
        verbose : ``boolean`` (default None)
            verbosity boolean flag. If left to None (default), use the reader verbosity flag.

        Returns
        -------
        dset : ``pymses.core.datasets.PointDataset``

        """
        amr_filename = filename_utils.amrlike_filename("amr", self._output_repos, self._iout, icpu, check_exists=True)

        # Verbosity
        if verbose is not None:
            if not isinstance(verbose, bool):
                raise AttributeError("'verbose' attribute must be a boolean value.")
            verb = verbose
        else:  # use reader verbosity flag
            verb = self._verbose

        # Load the AMR structure
        if verb:
            print("Reading amr data  : %s" % amr_filename)
        amr = read_ramses_amr_file(amr_filename, max_read_level=read_lmax, swap=self._swap)
        # amr_header, amr_struct = amr

        # Construct the octree
        dset = RamsesOctreeDataset(amr)
        for amr_file_type in list(self.ivars_descrs_by_file.keys()):
            if amr_file_type == "grav":
                is_hydro = 0
            else:
                is_hydro = 1
            ivars_to_read, descrs = self.ivars_descrs_by_file[amr_file_type]
            # Output file name
            fname = filename_utils.amrlike_filename(amr_file_type, self._output_repos, self._iout, icpu,
                                                    check_exists=False)

            # Read the data
            if verb:
                print("Reading %s : %s" % (amr_file_type.ljust(9), fname))
            data = read_ramses_amr_data_file(fname, is_hydro, amr, ivars_to_read, swap=self._swap,
                                             grav_compat=self._grav_compat)
            # try:
            #     data = read_ramses_amr_data_file(fname, is_hydro, amr, ivars_to_read, swap=self._swap)
            # except IOError as err:
            #     print err
            #     continue

            # Pack the fields into the octree thanks to the descriptors
            for desc in descrs:
                desc.gather(data, dset)

        # CPU id in 0-starting indexing
        dset.icpu = icpu - 1
        return dset


class CameraOctreeDataset(OctreeDataset):
    r"""
    Camera filtered octree dataset class

    contains all the relevant information about the AMR tree structure in a camera region.

    """

    def __init__(self, *args, **kwargs):
        if isinstance(args[0], tuple):  # Initialisation with amr dicts
            amr_dicts = args[0]
            super(CameraOctreeDataset, self).__init__(amr_dicts)
        else:
            assert isinstance(args[0], int) and isinstance(args[1], int)
            ngrid_max, ndim = args
            if "scalars" in kwargs:
                scalars = kwargs["scalars"]
            else:
                scalars = None
            if "vectors" in kwargs:
                vectors = kwargs["vectors"]
            else:
                vectors = None
            if "multivalued" in kwargs:
                mv = kwargs["multivalued"]  # must be a {'field_name': field_nvars, ...} dict !
            else:
                mv = None
            twotondim = (1 << ndim)
            twondim = 2 * ndim
            amr_head = {"ndim": ndim}
            amr_str = {"son_indices": -numpy.ones((ngrid_max, twotondim), dtype='i'),
                       "neighbors": -numpy.ones((ngrid_max, twondim), dtype='i'),
                       "grid_centers": numpy.zeros((ngrid_max, ndim)),
                       "cell_levels": numpy.zeros(ngrid_max, dtype='i'),
                       "ngrids": ngrid_max,
                       "o_igrid_max": 1}

            # Set center coords of the root grid
            amr_str["grid_centers"][0, :] = 0.5
            amr_str["cell_levels"][0] = 1

            super(CameraOctreeDataset, self).__init__((amr_head, amr_str))

            if scalars is not None:
                for fs in scalars:
                    self.add_scalars(fs, numpy.zeros((ngrid_max, twotondim)))
            if vectors is not None:
                for fv in vectors:
                    self.add_vectors(fv, numpy.zeros((ngrid_max, twotondim, ndim)))
            if mv is not None:
                for fmv_name, nvars in mv.items():
                    self.add_multivalued(fmv_name, numpy.zeros((ngrid_max, twotondim, nvars)))

        self.icpu = 0

    def restrict(self):
        """

        :return:
        """
        ng = self.amr_struct["o_igrid_max"]
        new_dset = CameraOctreeDataset(ng, self.amr_header["ndim"])
        new_dset.amr_struct["readlmax"] = self.amr_struct["readlmax"]
        new_dset.amr_struct["son_indices"] = self.amr_struct["son_indices"][:ng, :].copy()
        new_dset.amr_struct["neighbors"] = self.amr_struct["neighbors"][:ng, :].copy()
        new_dset.amr_struct["grid_centers"] = self.amr_struct["grid_centers"][:ng, :].copy()
        new_dset.amr_struct["cell_levels"] = self.amr_struct["cell_levels"][:ng].copy()
        for field in self.scalars:
            new_dset.add_scalars(field, self[field][:ng, :].copy())
        for field in self.vectors:
            new_dset.add_vectors(field, self[field][:ng, :, :].copy())
        for field in self.multivalued:
            new_dset.add_multivalued(field, self[field][:ng, :, :].copy())
        del self.amr_struct["o_igrid_max"]
        return new_dset


class CameraOctreeDatasource(RegionFilter):
    r"""Camera octree dataset

    This call the octree_build function which fills an octree structure with a RAMSES AMR dataset
    by taking into account the original AMR refinement only in the camera specific region. This
    regroup	cpu distributed octree datasets in a single octree dataset.

    Parameters
    ----------

    camera:
        camera that defines the region of interest for the octree dataset
    esize :
        extension size : extend the camera region
    ramses_amr_source
        ramses source to flatten
    radius (default None):
        define a sphere region like that : Sphere(camera.center, radius+esize)
        -> a RegionFilter is called on this region
    ngrid_max (default 2.000.000):
        Initial size of the created AMR array : 1e7 = 10 Millions <-> 3.8 GBytes of RAM memory
        This parameter has to be big enough to fit the local octree size.
    include_split_cells (default False):
        If True, the created octree will include all physical values of
        intermediary AMR resolution level (i.e. cells that are refined).
        If False, only leaf cell values are stored (this save memory
        and computation time for cell_to_points splatting rendering)
    """

    def __init__(self, camera, esize, ramses_amr_source, radius=None, ngrid_max=2000000,
                 include_split_cells=False):
        t0 = time()
        self.ngrid_max = int(ngrid_max)
        # Extended camera
        rs = numpy.ones((2, 3)) * esize
        rs[0, :] = -rs[0, :]
        from pymses.analysis.splatting.camera import ExtendedCamera

        ecam = ExtendedCamera(camera, rs)
        bb = ecam.get_bounding_box()

        if radius is None:
            # Init Filter with extended camera bounding box region
            reg = bb
        else:
            # Defined sphere region, checks the region includes the camera area
            zmax = numpy.max([camera.far_cut_depth, camera.distance])
            xmax = camera.region_size[0] / 2.
            ymax = camera.region_size[1] / 2.
            assert radius ** 2 >= (xmax ** 2 + ymax ** 2 + zmax ** 2)
            reg = Sphere(camera.center, radius + esize)

        super(CameraOctreeDatasource, self).__init__(reg, ramses_amr_source)

        # Max. read level setup
        lreq = camera.get_required_resolution()
        self.set_read_levelmax(lreq)

        # Init. self.dset to None
        self.dset = None
        self.build_dset(ecam, radius, include_split_cells=include_split_cells)
        print("CameraOctreeDatasource loaded up to level", lreq, \
            "with ngrids =", self.dset.amr_struct["ngrids"], \
            "(loading time = %.2fs" % (time() - t0), ")")

    def build_dset(self, camera, radius, include_split_cells=0):
        # Fetch camera info
        cc = camera.center
        ca = camera.get_camera_axis()
        corners = corner_points(numpy.zeros((1, 3)), numpy.ones(1))
        rcorners = camera.viewing_angle_rotation().transform_points(corners)
        bbmin = numpy.min(rcorners, axis=0)
        bbmax = numpy.max(rcorners, axis=0)
        cell_box = Box([bbmin, bbmax])
        cmb = camera.get_map_box()

        # Build octree dataset
        odset = None
        # Region-filtered octree structure creation
        for dset in self.iter_dsets():
            if odset is None:
                ndim = dset.amr_header["ndim"]
                multi_nvars = {fname: dset[fname].shape[-1] for fname in dset.multivalued}
                odset = CameraOctreeDataset(self.ngrid_max, ndim, scalars=dset.scalars, vectors=dset.vectors,
                                            multivalued=multi_nvars)
                odset.amr_struct["readlmax"] = dset.amr_struct["readlmax"]
            tree_utils.octree_build(odset, dset, (cc, ca, cmb, cell_box), radius,
                                    include_split_cells=include_split_cells)
            assert (odset.amr_struct["o_igrid_max"] < odset.amr_struct["ngrids"]), "Memory error : increase ngrid_max"

        # Neighbor indices computation
        tree_utils.octree_compute_neighbors(odset)

        # Set self.dset to new built octree
        # C'est moche, on restreint le tableau en memoire en redupliquant l'octree mais a la bonne taille
        self.dset = odset.restrict()

        # Set the data cpu list to [0]. The octree is built in memory, the source now only returns the dataset
        self.data_list = [0]

    def get_domain_dset(self, idomain, verbose=None):
        if self.dset is None:
            return RegionFilter.get_domain_dset(self, idomain, verbose)
        else:
            return self.dset


__all__ = ["RamsesOctreeDataset", "RamsesOctreeReader", "CameraOctreeDatasource", "CameraOctreeDataset"]
