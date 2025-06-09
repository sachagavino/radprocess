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

import wx
import os
import pymses


def open_about_box(frame):
    from wx.lib.wordwrap import wordwrap
    info = wx.AboutDialogInfo()
    info.SetName("About AMRViewer")
    info.SetVersion(pymses.__version__)
    info.SetCopyright("(C) CEA / IRFU 2009-2016")

    # Description
    # change the wx.ClientDC to use self.panel instead of self
    info.SetDescription(wordwrap(
            "The AMRViewer is a GUI based on PyMSES. It is designed to let RAMSES users who perform "
            "3D simulations have a quick and easy way to visualize the content of their simulation "
            "outputs. The visualization of a set of AMR physical/numerical quantities are enabled :\n\n"
            "Mass-weighted gas density maps\n"
            "Gas surface density maps\n"
            "Mass-weighted gas temperature maps\n"
            "Mass-weighted line-of-sight gas velocity maps\n"
            "Mass-weighted line-of-sight gas velocity dispersion maps\n\n"
            + pymses.__revision__,\
            350, wx.ClientDC(frame)))
    info.SetWebSite(("http://irfu.cea.fr/Projets/PYMSES", "PyMSES webpage"))
    info.AddDeveloper("Damien CHAPON")
    info.AddDeveloper("Marc LABADENS")

    # License text
    mp = os.path.dirname(globals()['__file__'])
    fname = os.path.join(mp, "LICENSE")
    file = open(fname)
    l = file.readlines(1000)
    licenseText = "".join(l)
    # change the wx.ClientDC to use self.panel instead of self
    info.SetLicense(wordwrap(licenseText, 500, wx.ClientDC(frame)))

    # Then we call wx.AboutBox giving it that info object
    wx.AboutBox(info)
