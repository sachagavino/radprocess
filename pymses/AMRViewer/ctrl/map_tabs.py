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

from ..view.widget_ids import WidgetIds as wid
from wx._core import PyDeadObjectError
from ..view.dialogs import *


class MapTabController(object):
    timing_ms = 20  # Refreshing time interval

    def __init__(self):
        self.active_tab = None
        self.tab_time = 0

    def NotebookChangeTab(self, event):
        """
        Changed notebook tab : update the name of the active tab
        """
        tab_id = event.GetSelection()
        if tab_id == 0:  # Region Finder tab
            active_tab = None
        else:
            WID = wid()
            active_tab = WID.TAB_NAME_LIST[tab_id - 1]
        self.ChangeActiveTab(active_tab)

    def ChangeActiveTab(self, active_tab):
        """
        Changed notebook tab : update the name of the active tab
        """
        self.active_tab = active_tab
        try:
            WID = wid()
            los = self.widgets[WID.EXPANDER_LOS]
            mag = self.widgets[WID.EXPANDER_GLASS]
            if self.active_tab is None:  # Region Finder tab
                # Expand los widget
                losCollapse = False
                self.SetRegionFinderHelpMessage()
            else:
                # Collapse los widget
                losCollapse = True
                self.gui.statusBar.SetStatusText(" Click and hold to measure distance between two points", 0)
            self.UpdateLogScale()  # update log_scale GUI view
            self.UpdateFraction()  # update fraction GUI view
            paramExpander = self.widgets[WID.EXPANDER_PARAMETERS]
            if self.active_tab == "velocity":
                paramExpander.cmap_cb.SetValue("RdBu_r")
                self.model.cmap = "RdBu_r"
            else:
                paramExpander.cmap_cb.SetValue("jet")
                self.model.cmap = "jet"
            los.Collapse(losCollapse)
            mag.Collapse(not losCollapse)
            if self.active_tab == "hdf5":
                self.gui.hdf5_cb.Show(True)
                self.gui.UpdateLayout(None)
            else:
                self.gui.hdf5_cb.Show(False)
            if self.active_tab == "stars_particles":
                self.gui.stars_age_instensity_dimmingCheckBox.Show(True)
                self.gui.UpdateLayout(None)
            else:
                self.gui.stars_age_instensity_dimmingCheckBox.Show(False)
            if self.active_tab == "levelmax" and self.model.refresh_levelmax_tab:
                tab = self.widgets[WID.TAB_IDS[self.active_tab]]
                tab.UpdateTabImage(verbose=True)
                tab.RefreshWholeTab(make_bitmap=True)
                self.model.refresh_levelmax_tab = False
            elif self.active_tab == "densityRT" and self.model.refresh_densityRT_tab:
                tab = self.widgets[WID.TAB_IDS[self.active_tab]]
                tab.UpdateTabImage(verbose=True)
                self.model.refresh_densityRT_tab = False
                tab.RefreshWholeTab(make_bitmap=True)

        except PyDeadObjectError as e:
            # avoid and error when leaving
            return

    def TabMapWindowOnMouse(self, event):
        """ Handles mouse events coming from a MapWindow of any MapTab
        of the notebook (right panel of the view)
        """
        WID = wid()
        win = event.GetEventObject()
        tab = self.widgets[WID.TAB_IDS[win.name]]

        # nothing to do if the window doesn't have a map to display
        if win.bitmap is None:
            return

        if event.Leaving():  # Leaving the window
            win.MoveCursor(None, None)
        else:
            t = event.GetTimestamp()
            if ((event.Moving() + event.Dragging()) * (
                        (t - self.tab_time) < MapTabController.timing_ms)):  # Refreshment frequency
                return
            self.tab_time = t

            # Cursor movement : Pointer position of the event
            i, j = event.GetPosition()
            win.MoveCursor(i, j)

            if event.LeftDown():
                # Left button pressed : SET behavior
                win.FreezeOrigin(i, j)

            if event.LeftUp():
                # Left button released : UNSET behavior
                win.FreezeOrigin(None, None)

        # Magnifier
        gle = self.widgets[WID.EXPANDER_GLASS]
        # Mouse wheel behavior to modify zoom box size
        wh = event.GetWheelRotation()
        # Wheel event to modify the Magnifier size
        if (wh != 0):
            wh = int(wh / abs(wh))
            index = gle.mag_factor_cb.GetSelection() + wh
            if index >= 0 and index < len(gle.mag_factor_list):
                gle.mag_factor_cb.Select(index)
                gle.set_mag_factor(index)

        # Display a magnified zone of the map in the Magnifier window
        gl_img = win.GetMagnifierBoxImage()
        gle.SetGlassImage(gl_img)

        # Display physical quantity value under the cursor
        pos = win.GetMapCursorPosition()
        self.model.SetGlassCenter(pos)
        val, unit = self.model.GetGlassValue(win.name)
        if ((val is None) + (unit is None)):
            label = "        ---        "
        else:
            if self.active_tab == "levelmax":
                label = "Level max = %i" % (round(val))
            else:
                label = "%.2e %s" % (val, unit)
        gle.RefreshPhysicalValue(label)
        # print the same information in the statusBar:
        self.gui.statusBar.SetStatusText(" Value : " + label, 1)
        rule_dist = win.GetRuleDist()
        if self.active_tab != "hdf5" and self.model.SetRuleDist(rule_dist):
            val, unit = self.model.GetRuleValue()
            if unit is not None:
                self.gui.statusBar.SetStatusText(" Measured distance : %.2f %s" % (val, unit), 1)
        tab.RefreshWholeTab()

    def TabMapWindowUpdateClick(self, event=None, do_all_processing=True, verbose=True):
        """
        Click on the update button of any TabMap
        """
        WID = wid()
        # Get the active tab
        tab = self.widgets[WID.TAB_IDS[self.active_tab]]
        if self.active_tab == "hdf5":
            if do_all_processing:
                win = HDF5FileDialog(self)
                win.run()
            if self.active_tab in self.model.image_computed and \
                            self.model.image_computed[self.active_tab] == True:
                tab.UpdateTabImage(verbose=verbose)
                tab.RefreshWholeTab(make_bitmap=True)
        elif self.model.ro is not None:
            # Update map
            if do_all_processing:
                nx, ny = tab.GetWindowSize()  # Size of the tab window
                self.model.UpdateTabMap(nx, self.active_tab)
            tab.UpdateTabImage(verbose=verbose)

            #			tab.MoveCursor(u=None, v=None, w=None)
            #			# Refresh cursor types on the tab window
            #			tab.RefreshWindowCursor()

            tab.RefreshWholeTab(make_bitmap=True)
