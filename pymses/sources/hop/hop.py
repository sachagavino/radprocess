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
from .file_formats import *
import os


class HOP(object):
    def __init__(self, input_file, hop_dir):
        r"""
        HOP simulation builder

        Parameters
        ----------

        input_file : ``string``

        hop_dir : ``string``
            HOP simulation directory path

        """
        assert os.path.isfile(input_file)
        self.input = input_file
        self.fname, ext = os.path.splitext(os.path.basename(self.input))
        assert os.path.isdir(hop_dir)
        self.dir = hop_dir
        self.root = os.path.join(self.dir, self.fname)
        self.hop_fname = "%s.hop" % self.root
        self.tag_fname = "%s.tag" % self.root
        self.den_fname = "%s.den" % self.root

    def process_hop(self, force=False):
        r"""

        Parameters
        ----------

        """
        if ((not force) * (os.path.isfile(self.hop_fname))):
            print("'%s' already exists." % self.hop_fname)
            return self.hop_fname
        if ((not force) * (os.path.isfile(self.den_fname))):
            print("'%s' already exists." % self.den_fname)
            hop_cmd = "hop -in %s -den %s -o %s" % (self.input,self.den_fname,self.root)
            os.system(hop_cmd)
            return self.den_fname
        hop_cmd = "hop -in %s -o %s" % (self.input, self.root)
        os.system(hop_cmd)
        return self.hop_fname

    def process_group(self, densities, force=False):
        r"""

        Parameters
        ----------

        """
        if ((not force) * (os.path.isfile(self.tag_fname))):
            print("'%s' already exists." % self.tag_fname)
            return self.tag_fname

        douter, dsaddle, dpeak = densities
        regroup_cmd = "regroup -root %s " % self.root + \
                      "-douter %f " % douter + \
                      "-dsaddle %f " % dsaddle + \
                      "-dpeak %f " % dpeak + \
                      "-o %s" % self.root
        os.system(regroup_cmd)
        return self.tag_fname

    def get_ngroups(self):
        r"""

        Parameters
        ----------


        """
        np, ng = HopTagHopFile(self.tag_fname).read_nparts_ngroups()
        return ng

    def iter_group(self, return_density=False):
        r"""

        Parameters
        ----------

        """
        tag_file = HopTagHopFile(self.tag_fname)
        group_ids = tag_file.read_group_id()

        if return_density:
            den_file = HopDenFile(self.den_fname)
            densities = den_file.read_den()
            for id_group in range(tag_file.ngroups):
                mask = (id_group == group_ids)
                dens = densities[mask]
                yield (id_group, mask, dens)
        else:
            for id_group in range(tag_file.ngroups):
                mask = (id_group == group_ids)
                yield (id_group, mask)

    def iter_group_hop(self, return_density=False):
        r"""

        Parameters
        ----------

        """
        hop_file = HopTagHopFile(self.hop_fname)
        group_ids = hop_file.read_group_id()

        if return_density:
            den_file = HopDenFile(self.den_fname)
            densities = den_file.read_den()
            for id_group in range(hop_file.ngroups):
                mask = (id_group == group_ids)
                dens = densities[mask]
                yield (id_group, mask, dens)
        else:
            for id_group in range(hop_file.ngroups):
                mask = (id_group == group_ids)
                yield (id_group, mask)

    def get_group_ids(self):
        r"""

        Parameters
        ----------

        """
        hop_file = HopTagHopFile(self.hop_fname)
        group_ids = hop_file.read_group_id()

        return group_ids


__all__ = ["HOP"]
