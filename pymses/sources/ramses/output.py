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
:mod:`pymses.sources.ramses.output` --- RAMSES output package
-------------------------------------------------------------

"""
from . import info, filename_utils, octree, particles, sources, sink
from .amr import read_ramses_amr_file
from .domain_decomposition import HilbertDomainDecomp
from .domain_decomposition import PlanarDomainDecomposition
from .field_descr import AmrDataFileFieldDescriptor
import pymses
from pymses.etc.rcsetup import rcConfiguration
from pymses.utils import constants as C

import os
import sys


class RamsesOutput(object):
    r"""
    Build a RamsesOutput object from an output repository path and the number of the output.
    Select first available output if no output numer is provided.

    Parameters
    ----------
    output_repos : ``string`` (default None)
        path to output repository
    iout : ``int``
        output number. If left to None (default), select the first available output detected in the directory.
    order : ``string`` of length 1 (default '=')
        byte order ('=' for native, '<' for little-endian, '>' for big-endian)
    verbose : ``boolean``
        verbosity boolean flag. Default True.
    check_endianness: ``bool``
        Performs data ENDIANNESS validity check. Default True.

    Examples
    --------

        >>> ro = RamsesOutput() # select current directory first output
        >>> ro = RamsesOutput("/data/Aquarius/outputs", 193)
    """
    def __init__(self, output_repos=None, iout=None, order='=', verbose=True, check_endianness=True):
        if output_repos is None:
            output_repos = "./"
        if iout is None:
            ilist = filename_utils.search_valid_outputs(output_repos)
            if len(ilist) == 0:
                raise AttributeError("'%s' is not a valid output repository.")
            output_repos = os.path.realpath(output_repos)
            iout = ilist[0]
        self._output_repos = output_repos
        self._iout = iout
        self._verbose = verbose
        self._local_amr_field_descr = None
        self._local_part_field_descr = None

        _output_path = self.output_path

        # Search for namelist files
        nmls = self.namelist_filenames
        if nmls is not None and pymses.hpc_config.is_root:
            print("  -> Found {nmls_number:d} namelist files (*.nml)".format(nmls_number=len(nmls)))

        # Check byte order
        if order not in ['=', '<', '>']:
            raise ValueError("Invalid byte order. It should be '=', '<' or '>'")
        self._swap = False
        if order != "=":
            sysorder = '<' if sys.byteorder == 'little' else '>'
            if order != sysorder:
                self._swap = True

        # Read the .info file
        self.info = info.read_ramses_info_file(self.info_filename)
        self._hydro_info = None
        self._ndim = self.info["ndim"]
        self._ncpu = self.info["ncpu"]

        if check_endianness:
            icpu = 1
            print("  -> Checking endianness for output #{ioutput:d}...".format(ioutput=self._iout))
            while True:
                try:
                    amr_filename = filename_utils.amrlike_filename("amr", self._output_repos, self._iout, icpu,
                                                                   check_exists=True)
                    try:
                        read_ramses_amr_file(amr_filename, swap=self._swap)
                        break
                    except IOError:
                        if self._swap:
                            new_endianness = sys.byteorder
                        else:
                            new_endianness = 'big' if sys.byteorder == 'little' else 'little'
                        print("     \_,-> Error, trying to switch to '{end:s}' endianness...".format(end=new_endianness))
                        self._swap = not self._swap
                        try:
                            read_ramses_amr_file(amr_filename, swap=self._swap)
                            break
                        except IOError as ioerr:
                            raise AttributeError("Data endiannes check failed. Try to change byte order. Error : "
                                                 "{err:s}".format(err=str(ioerr)))
                except IOError:  # AMR file not found for this CPU domain
                    if icpu == self._ncpu:
                        raise IOError("Cannot find any domain AMR file in this Ramses output directory !")
                    icpu += 1

        if not isinstance(verbose, bool):
            raise AttributeError("'verbose' attribute must be a boolean value.")

        self._cpu_list = list(range(1, self._ncpu + 1))
        if self.info["ordering"] == "hilbert":
            keys = self.info["dom_decomp_Hilbert_keys"]
            if self._verbose:
                print("  -> Computing hilbert minimal domain description for output #{ioutput:d}...".format(ioutput=iout))
            self.info["dom_decomp"] = HilbertDomainDecomp(self._ndim, keys[:-1], keys[1:],
                                                          (self.info["levelmin"], self.info["levelmax"]))
        elif self.info["ordering"] == 'planar':
            keys =self.info["dom_decomp_planar_keys"]
            self.info["dom_decomp"] = PlanarDomainDecomposition(self._ndim, keys)

        # List all the available files
        # self.output_files = filename_utils.output_files_dict(output_repos, iout)

        # Retrieve the domain decomposition
        self.dom_decomp = self.info["dom_decomp"]

        # Gather sink particle data, if any.
        try:
            self._sink_data_reader = sink.SinkParticleDataReader(self._output_repos, self._iout)
            if self._verbose:
                print("  -> Sink particle file found : '{sink_fname!s}'."\
                    .format(sink_fname=os.path.basename(self._sink_data_reader.sink_filepath)))
        except IOError:
            self._sink_data_reader = None

    @property
    def output_number(self):
        """
        Returns the simulation output number
        """
        return self._iout

    @property
    def output_repository(self):
        """
        Returns the path of the simulation output parent directory
        """
        return self._output_repos

    @property
    def output_path(self):
        """
        Returns the simulation output directory path
        """
        return filename_utils.output_path(self._output_repos, self._iout)

    @property
    def info_filename(self):
        return filename_utils.info_filename(self._output_repos, self._iout)

    @property
    def namelist_filenames(self):
        return filename_utils.namelist_filename_list(self._output_repos)

    @property
    def hydro_info(self):
        if self._hydro_info is None:
            hydro_fname = filename_utils.amrlike_filename("hydro", self._output_repos, self._iout, 1)
            self._hydro_info = info.read_ramses_hydro_header(hydro_fname, self._swap)
        return self._hydro_info

    def sink_data(self, sink_fields):
        return self._sink_data_reader.read(sink_fields)

    def box_size(self, unit=C.pc):
        """
        Returns the size of the simulation box in a user-defined unit.
        """
        if not isinstance(unit, C.Unit):
            raise AttributeError("'unit' must be a `pymses.utils.constants.unit.Unit` instance.")

        return self.info["unit_length"].express(unit)

    def _amr_fields(self):
        """
        Get the :class:`~pymses.sources.ramses.field_descr.AmrDataFileFieldDescriptor` instance, system default or
         user-defined instance.
        """
        if self._local_amr_field_descr is None:
            return rcConfiguration().Ramses.amr_fields
        else:
            return self._local_amr_field_descr

    def define_amr_scalar_field(self, *args, **kwargs):
        """
        Ignore the field description given by the pymsesrc files and add a new AMR scalar field to the
        :class:`~pymses.sources.ramses.output.RamsesOutput` field list.

        When adding an AMR field for the first time with one of the following method :

         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_scalar_field`,
         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_vector_field` or
         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_multivalued_field`

        an empty field list is initialised.

        Parameters
        ----------
        args : ``list``
            positional argument list forwarded to the
            :class:`~pymses.sources.ramses.field_descr.RamsesScalar` constructor
        kwargs: ``dict``
            keyword argument dictionary forwarded to the
            :class:`~pymses.sources.ramses.field_descr.RamsesScalar` constructor

        Examples
        --------

            >>> r = RamsesOutput("/data/simu/cosmo/run_fiducial/", 43)
            >>> r.define_amr_scalar_field("hydro", "rho", 0)
            >>> r.define_amr_scalar_field("hydro", "P", 4)
            >>> r.amr_fields()
            AmrDataFileFieldDescriptor : [RamsesScalar(name = "rho", ivar = 0, file = "hydro"),
            RamsesScalar(name = "P", ivar = 4, file = "hydro")]
        """
        if self._local_amr_field_descr is None:
            self._local_amr_field_descr = AmrDataFileFieldDescriptor()

        self._local_amr_field_descr.add_scalar(*args, **kwargs)

    def define_amr_vector_field(self, *args, **kwargs):
        """
        Ignore the field description given by the pymsesrc files and add a new AMR vector field to the
        :class:`~pymses.sources.ramses.output.RamsesOutput` field list.

        When adding an AMR field for the first time with one of the following method :

         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_scalar_field`,
         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_vector_field` or
         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_multivalued_field`

        an empty field list is initialised.

        Parameters
        ----------
        args : ``list``
            positional argument list forwarded to the
            :class:`~pymses.sources.ramses.field_descr.RamsesVector` constructor
        kwargs: ``dict``
            keyword argument dictionary forwarded to the
            :class:`~pymses.sources.ramses.field_descr.RamsesVector` constructor

        Examples
        --------

            >>> r = RamsesOutput("/data/simu/cosmo/run_fiducial/", 43)
            >>> r.define_amr_vector_field("hydro", "vel", [1, 2, 3])
            >>> r.define_amr_vector_field("grav", "g", [0, 1, 2])
            >>> r.amr_fields()
            AmrDataFileFieldDescriptor : [RamsesVector(name = "vel", ivars = [1, 2, 3], file = "hydro"),
            RamsesVector(name = "g", ivars = [0, 1, 2], file = "grav")]
        """
        if self._local_amr_field_descr is None:
            self._local_amr_field_descr = AmrDataFileFieldDescriptor()

        self._local_amr_field_descr.add_vector(*args, **kwargs)

    def define_amr_multivalued_field(self, *args, **kwargs):
        """
        Ignore the field description given by the pymsesrc files and add a new AMR multivalued field to the
        :class:`~pymses.sources.ramses.output.RamsesOutput` field list.

        When adding an AMR field for the first time with one of the following method :

         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_scalar_field`,
         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_vector_field` or
         * :meth:`~pymses.sources.ramses.output.RamsesOutput.define_amr_multivalued_field`

        an empty field list is initialised.

        Parameters
        ----------
        args : ``list``
            positional argument list forwarded to the :class:`~pymses.sources.ramses.field_descr.RamsesMultiValued`
            constructor
        kwargs: ``dict``
            keyword argument dictionary forwarded to the :class:`~pymses.sources.ramses.field_descr.RamsesMultiValued`
            constructor

        Examples
        --------

            >>> r = RamsesOutput("/data/simu/cosmo/run_fiducial/", 43)
            >>> r.define_amr_multivalued_field("hydro", "B", 4, 6)
            >>> r.amr_fields()
            AmrDataFileFieldDescriptor : [RamsesMultiValued(name = "B", ivars = [4, ..., 9], file = "hydro")]
        """
        if self._local_amr_field_descr is None:
            self._local_amr_field_descr = AmrDataFileFieldDescriptor()

        self._local_amr_field_descr.add_multivalued(*args, **kwargs)

    def amr_fields(self):
        """
        Print AMR field descriptor list
        """
        print(self._amr_fields())

    def amr_source(self, amr_read_fields, cpu_list=None, grav_compat=False, verbose=None):
        r"""
        Return a :class:`~pymses.sources.ramses.sources.RamsesAmrSource`, able to read a set of amr fields


        Parameters
        ----------
        amr_read_fields : ``list`` of ``strings``
            list of AMR data fields that needed to be read
        cpu_list : ``list`` of ``int`` (default None)
            If specified, restricts the cpu list to this list (for
            faster initialization). Default : cpu_list = range(1,ncpu+1)
        grav_compat; ``bool``
            Old Ramses versions (prior to commit 5dd90f3, 2012-10-04) and new Ramses versions (posterior to commit
            bce4454, 2015-07-09) work with nvar_file (number of scalar fields written in the file) parameter written in
            file header :
                \_,-> gravitational attraction vector field only (old version)
                \_,-> potential scalar field (phi) + gravitational attraction vector field (new version)

            For Ramses versions BETWEEN 2012-10-04 AND 2015-07-09, ndim + 1 variables were saved but only ndim was
            written	in the file header.
            To read output files written by this flavor of Ramses, set the `grav_compat` attribute to True.	Default
            False.
        verbose : ``boolean`` (default None)
            verbosity boolean flag. If left to None (default), use the RamsesOutput verbosity flag.

        Returns
        -------
        ramses_amr_source : :class:`~pymses.sources.ramses.sources.RamsesAmrSource`
            RAMSES AMR data source
        """
        # Init. reader
        ivars_descrs_by_file = self._amr_fields().gather_read_fields(amr_read_fields)
        if verbose is not None:
            if not isinstance(verbose, bool):
                raise AttributeError("'verbose' attribute must be a boolean value.")
            verb = verbose
        else:  # use RamsesOutput verbosity flag
            verb = self._verbose
        reader = octree.RamsesOctreeReader(self._output_repos, self._iout, ivars_descrs_by_file, self._swap,
                                           grav_compat, verbose=verb)

        if cpu_list is not None:
            data_list = cpu_list
        else:
            data_list = self._cpu_list

        return sources.RamsesAmrSource(reader, self.dom_decomp, data_list, self._ndim, amr_read_fields,
                                       (self.info["levelmin"], self.info["levelmax"]))

    def particle_source(self, field_list, cpu_list=None, long_int=False, select_stars=True, select_dark_matter=True,
                        verbose=None):
        r"""
        Return a RamsesParticleSource, able to read a set of user-defined particle data fields.


        Parameters
        ----------
        field_list : ``list`` of ``strings``
            list of particle data fields that need to be read
        cpu_list : ``list`` of ``int`` (default None)
            If specified, restricts the cpu list to this list (for
            faster initialization). Default : cpu_list = range(1,ncpu+1)
        long_int : ``bool``
            -DLONGINT Ramses compilation enabled ?
        select_stars :  ``boolean`` (default True)
            if True :  select and read STARS particles
                (with "epoch" field != 0)
        select_dark_matter :  ``boolean`` (default True)
            if True : select and read only DARK MATTER particles
                (with "epoch" field = 0)
        verbose : ``boolean`` (default None)
            verbosity boolean flag. If left to None (default), use the RamsesOutput verbosity flag.

        Returns
        -------
        ramses_part_source : :class:`~pymses.sources.ramses.sources.RamsesParticleSource`
            RAMSES particle data source

        """
        # Init. reader
        if verbose is not None:
            if not isinstance(verbose, bool):
                raise AttributeError("'verbose' attribute must be a boolean value.")
            verb = verbose
        else:  # use RamsesOutput verbosity flag
            verb = self._verbose
        reader = particles.RamsesParticleReader(self._output_repos, self._iout, field_list, self.info["boxlen"],
                                                select_stars=select_stars, select_dark_matter=select_dark_matter,
                                                swap=self._swap, long_int=long_int, verbose=verb)

        if cpu_list is not None:
            data_list = cpu_list
        else:
            data_list = self._cpu_list

        return sources.RamsesParticleSource(reader, self.dom_decomp, data_list, self._ndim, field_list,
                                            (self.info["levelmin"], self.info["levelmax"]))

    def cell_size(self, level="max", size_unit=None):
        """
        Returns the size of a cell of a given AMR level into a user-defined length unit, or in box unit if left to None.

        Parameters
        ----------
        level: ``string`` or  ``int``
            AMR level of refinement. 'min', 'max' or a positive integer value.
        size_unit: :class:`~pymses.utils.constants.unit.Unit` or None.
            unit of type 'length' or None (box unit). Default None

        Returns
        -------
        length: ``float``
            the conversion factor of the cell size expressed into the required size unit (or in box unit).

        Example
        -------
        >>> dx = ro.cell_size(12, size_unit=C.kpc)
        Level 12 : cell size = 0.24414 kpc
        >>> print dx
        0.24414
        """
        if level == "max":
            level = self.info["levelmax"]
        elif level == "min":
            level = self.info["levelmin"]
        else:
            if not isinstance(level, int) or level < 0:
                raise AttributeError("'level' must be 'min', 'max' or a positive integer value.")
            if level > self.info["levelmax"]:
                print("Warning : level attribute is greater than the max. AMR level of refinement. No cell will have " \
                      "that size")

        if size_unit is not None and not isinstance(size_unit, C.Unit):
            raise AttributeError("'size_unit' must be a `pymses.utils.constants.unit.Unit` instance or left to None.")

        # cell_size (in Box unit)
        cell_size = 1. / 2**level

        if size_unit is not None:
            cell_size *= self.info["unit_length"].express(size_unit)
        return cell_size

    def show_cell_size(self, level, size_unit=None):
        """
        Print the size of a cell of a given AMR level.

        Parameters
        ----------
        level: ``int``
            AMR level of refinement
        size_unit: :class:`~pymses.utils.constants.unit.Unit`
            unit of type 'length'

        See also
        --------
        :meth:`~pymses.sources.ramses.output.RamsesOutput.cell_size()` method.
        """
        length = self.cell_size(level, size_unit=size_unit)
        if size_unit is None:
            print("Level %i : cell size = %g (in box unit)" % (level, length))
        else:
            print("Level %i : cell size = %g %s" % (level, length, size_unit.name))


__all__ = ["RamsesOutput"]
