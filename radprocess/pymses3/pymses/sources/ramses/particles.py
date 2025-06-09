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
:mod:`pymses.sources.ramses.particles` --- RAMSES output particle file reading module
-------------------------------------------------------------------------------------
"""
import os
import numpy
from . import _read_ramses
from . import filename_utils
from pymses.core.datasets import PointDataset
from pymses.core.reader import MultiFileDataReader


def read_ramses_particle_file(part_fname, lbox, swap=False, long_int=False):
    """
    Reads a RAMSES .part file into memory

    Parameters
    ----------

    part_fname : ``string``
        filename of the .part file
    lbox : ``float``
        simulation box width
    swap : ``bool``
        swap bytes if True (to match different endianness)
    long_int  : ``bool``
        -DLONGINT Ramses compilation option enabled ?

    Returns
    -------
    header_dict : ``dict``
        particle file header dictionary
    part_dict : ``dict``
        particle data dictionary
    """

    try:
        header, part = _read_ramses.read_parts(part_fname, lbox, int(swap), int(long_int))
    except _read_ramses.RamsesIOError as err:
        raise IOError(str(err))

    # Extract icpu from filename
    try:
        icpu = int(os.path.basename(part_fname).split("out")[1])
    except ValueError:
        raise RuntimeError("cannot parse icpu from file name")
    header["icpu"] = icpu

    return header, part


class RamsesParticleReader(MultiFileDataReader):
    """
    RAMSES particle file reader

    Parameters
    ----------
    output_repos              : :  ``string``
        path to output repository
    iout          : : ``int``
        output number
    field_list              : :  ``list of string``
        fields to read like ["mass", "level", "epoch"]
    lbox          : : ``float``
        box length (particules positions are divided by that value)
    fortran_file_kwargs :  ``dictionary`` (default None)
        useless here
    select_stars :  ``boolean`` (default True)
        if True :  select and read STARS particules
            (with "epoch" field != 0)
    select_dark_matter :  ``boolean`` (default True)
        if True : select and read only DARK MATTER particules
            (with "epoch" field = 0)
    swap : ``boolean`` (default False)
        if True : swap bytes (to match different endianness)
    long_int : ``boolean`` (default False)
        -DLONGINT Ramses compilation option enabled ?
    verbose :  ``boolean`` (default True)
        print a line in console for every file read
    """

    def __init__(self, output_repos, iout, field_list, lbox, fortran_file_kwargs=None, select_stars=True,
                 select_dark_matter=True, swap=False, long_int=False, verbose=True):
        self._output_repos = output_repos
        self._iout = iout
        self.field_list = field_list
        self.lbox = lbox

        if fortran_file_kwargs is None:
            self._ffkwargs = {}
        else:
            self._ffkwargs = fortran_file_kwargs
        self.select_stars = select_stars
        self.select_dark_matter = select_dark_matter
        self.swap = swap
        self._long_int = long_int
        self._verbose = verbose

    def read_file(self, icpu, read_lmax=None, verbose=None):
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
        part_fname = filename_utils.amrlike_filename("part", self._output_repos, self._iout, icpu, check_exists=True)
        keep_epoch_field_after_filter = True

        if (not self.select_stars or not self.select_dark_matter) and ("epoch" not in self.field_list):
            keep_epoch_field_after_filter = False
            self.field_list.append("epoch")

        # Read the particles data from the particle file
        if verbose is not None:
            if not isinstance(verbose, bool):
                raise AttributeError("'verbose' attribute must be a boolean value.")
            verb = verbose
        else:  # use reader verbosity flag
            verb = self._verbose
        if verb:
            print("Reading particles : %s" % part_fname)
        part_header, part_data = read_ramses_particle_file(part_fname, self.lbox, self.swap, self._long_int)

        if "epoch" not in part_data:
            if "epoch" in self.field_list:
                self.field_list.remove("epoch")
            if not self.select_stars or not self.select_dark_matter:
                print("Warning : \"epoch\" field not found -> no particle selection is done here")
                self.select_stars = True
                self.select_dark_matter = True
        if not self.select_stars or not self.select_dark_matter:
            # We have to filter stars or dm parts
            mask = numpy.ones(part_data["pos"].shape[0], 'bool')
            if not self.select_stars:
                mask *= part_data["epoch"] == 0  # select dm only
            if not self.select_dark_matter:
                mask *= part_data["epoch"] != 0  # select stars only
            if part_data["pos"].shape[0] != 0:
                dset = PointDataset(part_data["pos"][mask])  # Dataset creation
                # Filling Dataset with the user-chosen fields
                if "vel" in self.field_list:
                    dset.add_vectors("vel", part_data["vel"][mask])
                if "mass" in self.field_list:
                    dset.add_scalars("mass", part_data["mass"][mask])
                if "id" in self.field_list:
                    dset.add_scalars("id", part_data["id"][mask])
                if "level" in self.field_list:
                    dset.add_scalars("level", part_data["level"][mask])
                if "epoch" in self.field_list and keep_epoch_field_after_filter:
                    dset.add_scalars("epoch", part_data["epoch"][mask])
                if "metal" in self.field_list:
                    dset.add_scalars("metal", part_data["metal"][mask])
                return dset

        # Else, we don't have to filter stars or dm parts
        dset = PointDataset(part_data["pos"])  # Dataset creation
        # Filling Dataset with the user-chosen fields
        if "vel" in self.field_list:
            dset.add_vectors("vel", part_data["vel"])
        if "mass" in self.field_list:
            dset.add_scalars("mass", part_data["mass"])
        if "id" in self.field_list:
            dset.add_scalars("id", part_data["id"])
        if "level" in self.field_list:
            dset.add_scalars("level", part_data["level"])
        if "epoch" in self.field_list and keep_epoch_field_after_filter:
            dset.add_scalars("epoch", part_data["epoch"])
        if "metal" in self.field_list:
            dset.add_scalars("metal", part_data["metal"])

        return dset
