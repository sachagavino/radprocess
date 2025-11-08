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


class WidgetIds(object):
    class __wid(object):
        def __init__(self):
            self.ANY = wx.ID_ANY

            # Tools
            self.TOOL_OPEN = wx.NewId()
            self.TOOL_RESET = wx.NewId()
            self.TOOL_PREVIOUS = wx.NewId()
            self.TOOL_QUIT = wx.NewId()
            self.TOOL_UPDATE = wx.NewId()
            self.TOOL_BUILD_SOURCE = wx.NewId()
            self.TOOL_SAVE_IMAGE = wx.NewId()
            self.TOOL_GLNEMO2 = wx.NewId()
            self.TOOL_CLOSE = wx.NewId()

            # Menu
            self.MENU_OPEN = wx.ID_OPEN
            self.MENU_SAVECAM = wx.NewId()
            self.MENU_LOADCAM = wx.NewId()
            self.MENU_QUIT = wx.ID_EXIT
            self.MENU_ABOUT = wx.ID_ABOUT

            # Expanders
            self.EXPANDER_RAMSES = wx.NewId()
            self.EXPANDER_PARAMETERS = wx.NewId()
            self.EXPANDER_LOS = wx.NewId()
            self.EXPANDER_GLASS = wx.NewId()

            # Line-of-sight widgets
            self.LOS_A3D = wx.NewId()
            self.LOS_TOOL_X = wx.NewId()
            self.LOS_TOOL_MX = wx.NewId()
            self.LOS_TOOL_Y = wx.NewId()
            self.LOS_TOOL_MY = wx.NewId()
            self.LOS_TOOL_Z = wx.NewId()
            self.LOS_TOOL_MZ = wx.NewId()
            self.LOS_TOOL_ROT_LEFT = wx.NewId()
            self.LOS_TOOL_ROT_RIGHT = wx.NewId()
            self.LOS_TOOL_ROT_TOP = wx.NewId()
            self.LOS_TOOL_ROT_BOTTOM = wx.NewId()

            # Magnifying glass widgets
            self.GLASS_WIN = wx.NewId()
            self.GLASS_VAL = wx.NewId()

            # Region finder tab
            self.RF_TAB = wx.NewId()
            self.RF_TAB_LOS_WIN = wx.NewId()
            self.RF_TAB_U_WIN = wx.NewId()
            self.RF_TAB_V_WIN = wx.NewId()

            # Other tabs
            # Before any addition to these 2 lists, let's say "item", make sure :
            # - you added a xx*20 pixels image in the view/icons directory named 'item.png'
            self.TAB_NAME_LIST = ["levelmax", "densityRT", "density", "Sigma",  # "dm_particles", "stars_particles",
                                  "temperature", "velocity", "custom", "hdf5"]  # "transfer_function"
            self.TAB_DESCR_LIST = {
                "levelmax": "Max. AMR level of refinement along the line-of-sight",
                "densityRT": "Density ray tracing",
                "density": "Mass-weighted gas density",
                "Sigma": "Gas surface density",
                "temperature": "Mass-weighted gas temperature",
                "velocity": "Mass-weighted line-of-sight gas velocity",
                "dm_particles": "Dark matter particles mass, with size of splatting points corresponding to particle level",
                "stars_particles": "Stars particles mass, with size of splatting points corresponding to particle level",
                "transfer_function": "Transfer function ray tracing example",
                "custom": "User defined custom map",
                "hdf5": "PYMSES image plot utils HDF5 file viewer"}
            ####### DO NOT MODIFY BELOW THIS LINE ##########
            self.TAB_WIN_IDS = {}
            self.TAB_IDS = {}
            self.TAB_MENU_IDS = {}
            i = 1
            self.Tab_Number_for_selection = {}
            for tab_name in self.TAB_NAME_LIST:
                self.TAB_WIN_IDS[tab_name] = wx.NewId()
                self.TAB_IDS[tab_name] = wx.NewId()
                self.TAB_MENU_IDS[tab_name] = wx.NewId()
                self.Tab_Number_for_selection[tab_name] = i
                i += 1
                ################################################

    instance = None

    def __new__(c):  # _new_ is always a class method
        if not WidgetIds.instance:
            WidgetIds.instance = WidgetIds.__wid()
        return WidgetIds.instance

    def __getattr__(self, attr):
        return getattr(self.instance, attr)

    def __setattr__(self, attr, val):
        pass
