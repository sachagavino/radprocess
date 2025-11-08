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
from .widget_ids import WidgetIds as wid
from .region_finder import RegionFinderTab
from .map_tabs import MapTab
from .utils import get_icon
from sys import platform as _platform
import os


class Notebook(wx.Notebook):
    def __init__(self, parent, model, controller):
        WID = wid()
        super(Notebook, self).__init__(parent, WID.ANY)
        self.rf_tab = RegionFinderTab(self, model, controller)

        # Physical/numerical quantities map tabs
        self.map_tabs = dict.fromkeys(WID.TAB_NAME_LIST)

        # OS specific test
        if _platform == "darwin":
            is_mac_OS = True  # OS X user : specific icons
        else:
            is_mac_OS = False  # _platform may be "linux" "linux2" or "win32"

        il = wx.ImageList(25, 25)
        icon_id = {}
        for name in list(self.map_tabs.keys()):
            self.map_tabs[name] = MapTab(self, name, model, controller)
            fname_mac = get_icon("%s_mac.png" % name)
            if is_mac_OS and os.path.isfile(fname_mac):
                b = wx.Bitmap(fname_mac)
            else:
                b = wx.Bitmap(get_icon("%s.png" % name))
            icon_id[name] = il.Add(b)
        self.AssignImageList(il)

        self.model = model

        self.AddPage(self.rf_tab, "Region finder")
        for index in range(len(WID.TAB_NAME_LIST)):
            name = WID.TAB_NAME_LIST[index]
            self.AddPage(self.map_tabs[name], "")
            self.map_tabs[name].Hide()
            self.SetPageImage(index + 1, icon_id[name])

        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, controller.NotebookChangeTab)

    def reset(self):
        self.rf_tab.reset()
        for tab in list(self.map_tabs.values()):
            tab.reset()
