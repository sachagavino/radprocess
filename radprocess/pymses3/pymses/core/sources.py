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
"""
:mod:`pymses.core.sources` --- PyMSES generic data source module
================================================================

"""
from time import time


class DataSource(object):
    r"""
    Base class for all data source objects

    Parameters
    ----------
    dom_decomp : ````
        domain decomposition object instance
    data_list : ``list``
        list of the data domain indices
    field_list : ``list``
        list of the field names provided by the datasource
    read_lmax : ``int `
        higher resolution level of the data read by the datasource
    ndim : ``int``
        number of dimensions of the data provided by the datasource

    """
    AMR_SOURCE = 10
    PARTICLE_SOURCE = 11

    def __init__(self, dom_decomp=None, data_list=None, field_list=None, read_lmax=None, ndim=None):
        self.dom_decomp = dom_decomp
        if data_list is None:
            self.data_list = []
        else:
            if not isinstance(data_list, list):
                raise AttributeError("'data_list' must be a list object. Got '%s'." % type(data_list))
            self.data_list = data_list

        if field_list is None:
            self.field_list = []
        else:
            if not isinstance(field_list, list):
                raise AttributeError("'field_list' must be a list object. Got '%s'." % type(field_list))
            self.field_list = field_list

        self._read_lmax = read_lmax
        self.ndim = ndim

    def __del__(self):
        del self.dom_decomp
        del self.data_list

    @property
    def ndataset(self):
        return len(self.data_list)

    def iter_idata(self):
        for idata in self.data_list:
            yield idata

    def get_read_levelmax(self):
        """
        Returns
        -------
        Get the maximum resolution level value of the data read by the datasource
        """
        return self._read_lmax

    def set_read_levelmax(self, read_lmax):
        r"""Sets the maximum resolution level to read in the datasource.

        Parameters
        ----------
        read_lmax : ``int``
            max. resolution level to read

        """
        lmax = int(read_lmax)
        if lmax <= 0:
            raise AttributeError("'read_lmax' parameter value is not valid. A positive integer value is required.")
        self._read_lmax = lmax

    def source_type(self):
        raise NotImplementedError()

    def iter_dsets(self, verbose=None):
        r"""
        Datasets iterator method. Yield datasets from the datasource

        Parameters
        ----------
        verbose: ``bool``
            verbosity boolean flag. Default None.
        """
        for idata in self.data_list:
            yield self.get_domain_dset(idata, verbose)

    def __iter__(self):
        return self.iter_dsets()

    def get_domain_dset(self, idomain, verbose=None):
        raise NotImplementedError()

    def flatten(self, verbose=None):
        """
        Read each data file and concatenate resulting dsets. This method tries to use multiprocessing if possible.

        Returns
        -------
        fdset : flattened dataset
        """
        startingTime = time()
        dsets = list(self.iter_dsets(verbose=verbose))
        # TODO : proper multiprocessing management
        print("Read and filter time : %.2f s" % (time() - startingTime))

        if len(dsets) == 0:
            return None
        else:
            return dsets[0].concatenate(dsets)


class Filter(DataSource):
    r"""
    Data filter generic class.

    Parameters
    ----------
    source: base data source
    """
    def __init__(self, source):
        super(Filter, self).__init__(source.dom_decomp, source.data_list, source.field_list, ndim=source.ndim)
        self._source = source

    @property
    def source(self):
        return self._source

    def __del__(self):
        super(Filter, self).__del__()
        del self._source

    def set_read_levelmax(self, max_read_level):
        r"""
        [Overwritten] Apply the set_read_levelmax() method to the filter's source.

        Parameters
        ----------
        max_read_level : ``int``
            max. resolution level of the data provided by the filter

        """
        self._source.set_read_levelmax(max_read_level)

    def get_read_levelmax(self):
        """
        Returns
        -------
        Get the maximum resolution level value of the data read by the filter's base datasource
        """
        return self._source.get_read_levelmax()

    def source_type(self):
        r"""
        Returns
        -------
        type : ``int``
            the result of the : class:`~pymses.cors.sources.Datasource.source_type()` method of the filter's source.

        """
        return self._source.source_type()

    def filtered_dset(self, dset):
        r"""
        Abstract `filtered_dset()` method

        """
        raise NotImplementedError()

    def get_domain_dset(self, idomain, verbose=None):
        r"""
        Get the filtered result of `self.source.get_domain_dset(idomain)`

        Parameters
        ----------
        idomain : ``int``
            number of the domain from which the data is required

        Returns
        -------
        dset : ``AbstractDataset``
            the filtered dataset corresponding to the given `idomain`

        """
        return self.filtered_dset(self._source.get_domain_dset(idomain, verbose=verbose))


class SubsetFilter(Filter):
    r"""
    SubsetFilter class. Selects a subset of datasets to read from the
    datasource

    Parameters
    ----------
    data_sublist : ``list`` of ``int``
        list of the selected dataset index to read from the datasource

    """

    def __init__(self, data_sublist, source):
        super(SubsetFilter, self).__init__(source)
        if data_sublist is not None:
            self.data_list = [i for i in data_sublist if i in self.data_list]

    def filtered_dset(self, dset):
        return dset


__all__ = ["DataSource", "Filter", "SubsetFilter"]
