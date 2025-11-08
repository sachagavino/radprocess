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
:mod:`pymses.sources.ramses.sources` --- RAMSES data sources module
-------------------------------------------------------------------

"""
from pymses.core import DataSource


class RamsesGenericSource(DataSource):
    r"""
    RAMSES generic data source

    Parameters
    ----------
    reader: TODO
    dom_decomp: TODO
    cpu_list:``list`` of ``int``
        list of cpu index
    ndim: ``int``
        number of dimensions
    field_list: ``list`` of ``string``
        list of field names to read
    data_amr_levels: (levelmin, levelmax) ``tuple``
        tuple containing the min. and max. AMR data structure levels of refinement
    """
    def __init__(self, reader, dom_decomp, cpu_list, ndim, field_list, data_amr_levels=None):
        if data_amr_levels is None:
            self._amrdata_levelmin = 1
            self._amrdata_levelmax = 999
        else:
            if not isinstance(data_amr_levels, tuple) and not isinstance(data_amr_levels, list):
                raise AttributeError("'data_amr_elevels' attribute is not a valid (levelmin, levelmax) tuple value.")
            if not isinstance(data_amr_levels[0], int) or not isinstance(data_amr_levels[1], int):
                raise AttributeError("'data_amr_elevels' attribute must contain integer values : (levelmin, levelmax).")

            self._amrdata_levelmin, self._amrdata_levelmax = data_amr_levels

        # Build datasource with read_lmax = max. AMR level available in the data
        super(RamsesGenericSource, self).__init__(dom_decomp, cpu_list, field_list, ndim=ndim,
                                                  read_lmax=self._amrdata_levelmax)
        self.reader = reader

    def get_domain_dset(self, icpu, verbose=None):
        r"""
        Data source reading method

        Parameters
        ----------
        icpu : ``int``
            CPU file number to read
        verbose: ``bool`` or None.
            verbosity boolean flag. Default None.

        Returns
        -------
        dset : ``AbstractDataset``
            the dataset containing the data from the given cpu number file

        """
        return self.reader.read_file(icpu, self.get_read_levelmax(), verbose=verbose)

    @property
    def amr_data_levelmax(self):
        return self._amrdata_levelmax

    @property
    def amr_data_levelmin(self):
        return self._amrdata_levelmin


class RamsesAmrSource(RamsesGenericSource):
    r"""
    RAMSES AMR data source class

    """

    def __init__(self, reader, dom_decomp=None, cpu_list=None, ndim=None, field_list=None, data_amr_levels=None):
        super(RamsesAmrSource, self).__init__(reader, dom_decomp, cpu_list, ndim, field_list, data_amr_levels)

    def source_type(self):
        r"""
        Returns
        -------
        t: DataSource.AMR_SOURCE

        """
        return DataSource.AMR_SOURCE

    def set_read_levelmax(self, read_lmax):
        r"""Sets the maximum resolution level to read in the datasource.

        Parameters
        ----------
        read_lmax : ``int``
            max. resolution level to read

        """
        lmax = int(read_lmax)

        # Want to set read_levelmax to an AMR level greater than max. available AMR level in the data => no effect
        if lmax > self._amrdata_levelmax:
            return

        # Want to set read_levelmax to an AMR level lower (strictly) than max. available AMR  level in the data while
        # coarse grid values are not valid => no effect. Force reading AMR data down to max. available level !
        if lmax < self._amrdata_levelmax and not self.valid_coarse_grid_values:
            return

        # Call base method
        super(RamsesAmrSource, self).set_read_levelmax(read_lmax)

    @property
    def valid_coarse_grid_values(self):
        """
        Returns
        -------
        b: ``bool``
            True when levelmin < levelmax. WARNING : when levelmin==levelmax, AMR field values for ilevel<levelmin
            corresponds to initial conditions and are not updated afterwards by Ramses !!!
            Otherwise return False (levelmax==levelmin).
        """
        valid_coarse_grid_values = self.amr_data_levelmin < self.amr_data_levelmax
        return valid_coarse_grid_values


class RamsesParticleSource(RamsesGenericSource):
    r"""
    RAMSES particle data source class

    """

    def __init__(self, reader, dom_decomp=None, cpu_list=None, ndim=None, field_list=None, data_amr_levels=None):
        super(RamsesParticleSource, self).__init__(reader, dom_decomp, cpu_list, ndim, field_list, data_amr_levels)

    def source_type(self):
        r"""
        Returns
        -------
        DataSource.PARTICLE_SOURCE

        """
        return DataSource.PARTICLE_SOURCE


__all__ = ["RamsesGenericSource", "RamsesAmrSource", "RamsesParticleSource"]
