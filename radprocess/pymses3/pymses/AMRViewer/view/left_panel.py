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
from wx.lib.scrolledpanel import ScrolledPanel
from .expander import *


class LeftPanel(ScrolledPanel):
    def __init__(self, parent, model, controller):
        super(LeftPanel, self).__init__(parent, wx.ID_ANY, size=(245, -1))

        # Sizer
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(sizer)

        # Parameters expander
        self.cpp = ParametersExpander(self, model, controller)
        sizer.Add(self.cpp, 0, wx.ALL | wx.EXPAND, 1)

        # Line-of-sight axis expander
        self.cpa = LineOfSightExpander(self, model, controller)
        sizer.Add(self.cpa, 0, wx.ALL | wx.EXPAND, 1)

        # Magnifying glass expander
        self.cpg = MagGlassExpander(self, model, controller)
        sizer.Add(self.cpg, 0, wx.ALL | wx.EXPAND, 1)

        self.SetupScrolling(scroll_x=False)

        # Left Panel
        size = self.GetSize()
        self.SetMinSize((size[0], -1))
        self.SetMaxSize((size[0], -1))

    def reset(self):
        self.cpp.reset()
        self.cpa.reset()
        self.cpg.reset()
