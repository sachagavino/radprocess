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
from pymses.core.datasets import PointDataset
import os


class HopDenFile(object):
    def __init__(self, fname):
        assert os.path.isfile(fname), "HOP .den file '%s' doesn't exist." % fname
        self.fname = fname

    def read_den(self):
        self._ff = file(self.fname, 'r')
        self.nden = numpy.fromfile(self._ff, count=1, dtype='i')
        densities = numpy.fromfile(self._ff, count=self.nden, dtype='f')
        self._ff.close()
        return densities


class HopTagHopFile(object):
    def __init__(self, fname, **kwargs):
        assert os.path.isfile(fname), "HOP .hop/.tag file '%s' doesn't exist." % fname
        self.fname = fname

    def read_nparts_ngroups(self, close_after=True):
        self._ff = file(self.fname, 'r')
        self.nparts = numpy.fromfile(self._ff, count=1, dtype='i')
        self.ngroups = numpy.fromfile(self._ff, count=1, dtype='i')
        if close_after:
            self._ff.close()
        return (self.nparts, self.ngroups)

    def read_group_id(self):
        np, ng = self.read_nparts_ngroups(close_after=False)
        group_ids = numpy.fromfile(self._ff, count=self.nparts, dtype='i')
        self._ff.close()
        return group_ids


class HopPosFile(object):
    def __init__(self, fname, **kwargs):
        assert os.path.isfile(fname), "HOP pos. file '%s' doesn't exist." % fname
        self.fname = fname

    def read(self):
        dat = numpy.loadtxt(self.fname)
        ds = PointDataset(dat[:, 3:6])
        ds.add_scalars("rank", numpy.array(dat[:, 0], dtype='i'))
        ds.add_scalars("npart", numpy.array(dat[:, 1], dtype='i'))
        ds.add_scalars("m", dat[:, 2])
        ds.add_vectors("vel", dat[:, 6:9])
        return ds


__all__ = ["HopDenFile", "HopTagHopFile", "HopPosFile"]
