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

from .ramses_sim import RamsesSimController
from .los_axes3d import LineOfSightController
from .region_finder import RegionFinderController
from .map_tabs import MapTabController
from ..model.ramses import map_engine
from ..view.dialogs import SelectParametersDialog
from ..view.widget_ids import WidgetIds as wid
import wx
from ..model.model import Model
from threading import *
import csv, os
from pymses.utils import misc


class Controller(RamsesSimController, LineOfSightController, RegionFinderController, MapTabController):
    def __init__(self, gui):
        self.gui = gui
        RamsesSimController.__init__(self)
        LineOfSightController.__init__(self)
        RegionFinderController.__init__(self)
        MapTabController.__init__(self)
        self.model = Model()
        self.log_scale_dict = {}
        self.fraction_dict = {}
        self.default_fraction = "0.1"

    def SetWidgetsDict(self, widgets):
        self.widgets = widgets

    def OnReset(self, event):
        self.model.reset()
        self.gui.reset()
        self.log_scale_dict = {}
        self.fraction_dict = {}
        self.LoadConfigFromCSV()

    def OnPrevious(self, event):
        WID = wid()
        if self.active_tab is not None:
            self.ChangeActiveTab(None)
            self.gui.rp.ChangeSelection(0)
        if len(self.model.map_name_list_history) > 1:
            self.model.region_finder_previous_cam()
            self.gui.reset()
            self.UpdateImages(do_all_processing=False)

    def SaveLosCameraToFile(self, h5fname):
        return self.model.SaveLosCameraToFile(h5fname)

    def LoadLosCameraFromFile(self, h5fname):
        self.model.reset()
        self.model.LoadLosCameraFromFile(h5fname)
        self.gui.reset()

    def LoadmapFromHDF5(self, h5fname):
        self.model.LoadmapFromHDF5(h5fname)

    def SaveConfigInCSV(self):
        """
        Saves the config parameters into a csv (Comma Separated Values) file.
        Use ~/.pymses.cfg
        """
        try:
            paramExpander = self.widgets[wid().EXPANDER_PARAMETERS]
            config_file_path = os.path.expanduser("~") + "/.pymses.cfg"
            # print "Saving default parameters in ", config_file_path
            with open(config_file_path, 'wb') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=',')
                spamwriter.writerow(['PyMSES parameter csv file (this ordering is required)'])
                spamwriter.writerow(['Region finder map', paramExpander.regionFinderMapType_cb.GetValue()])
                spamwriter.writerow(['NUMBER_OF_PROCESSES_LIMIT', paramExpander.Multiprocessing_cb.GetValue()])
                spamwriter.writerow(['rf multiprocessing', paramExpander.MultiprocessingBtn.GetValue()])
                spamwriter.writerow(['resolution', paramExpander.lower_resolution_cb.GetValue()])
                spamwriter.writerow(['fft_factor', paramExpander.FFTkernelSize_cb.GetValue()])
                spamwriter.writerow(['remember_data', paramExpander.rememberSomeDataBtn.GetValue()])
                spamwriter.writerow(['random_shift', paramExpander.random_shiftCheckBox.GetValue()])
                spamwriter.writerow(['colormap', paramExpander.cmap_cb.GetValue()])
                spamwriter.writerow(['gaussian_blur', paramExpander.AdaptiveGaussianBlurBtn.GetValue()])
                spamwriter.writerow(['log_scale', paramExpander.log_CheckBox.GetValue()])
                spamwriter.writerow(['fraction', paramExpander.fraction_cb.GetValue()])
                spamwriter.writerow(['GUI_thread', paramExpander.ThreadCheckBox.GetValue()])

        except:
            None

    def LoadConfigFromCSV(self):
        """
        Load the config parameters from a csv (Comma Separated Values) file.
        Use ~/.pymses.cfg
        """
        try:
            def bool_from_str(str):
                if str != "False":
                    return True
                else:
                    return False

            paramExpander = self.widgets[wid().EXPANDER_PARAMETERS]
            config_file_path = os.path.expanduser("~") + "/.pymses.cfg"
            print("Loading default parameters from ", config_file_path)
            with open(config_file_path, 'rb') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                next(spamreader), " intro comment ! "
                rfMapType = spamreader.next()[1]
                paramExpander.regionFinderMapType_cb.SetValue(rfMapType)
                self.model.regionFinderMapType = rfMapType
                process_limit = spamreader.next()[1]
                paramExpander.Multiprocessing_cb.SetValue(process_limit)
                # misc.NUMBER_OF_PROCESSES_LIMIT = int(process_limit)
                multiprocessing = bool_from_str(spamreader.next()[1])
                paramExpander.MultiprocessingBtn.SetValue(multiprocessing)
                self.model.multiprocessing = multiprocessing
                resolution = spamreader.next()[1]
                paramExpander.lower_resolution_cb.SetValue(resolution)
                paramExpander.update_lower_resolution_cb(resolution=resolution, save_CSV=False)
                size = spamreader.next()[1]
                paramExpander.FFTkernelSize_cb.SetValue(size)
                paramExpander.updateFFTkernelSize(size=size, save_CSV=False)
                remember = bool_from_str(spamreader.next()[1])
                paramExpander.rememberSomeDataBtn.SetValue(remember)
                self.model.rememberSomeData = remember
                rdm = bool_from_str(spamreader.next()[1])
                paramExpander.random_shiftCheckBox.SetValue(rdm)
                self.model.random_shift = rdm
                cmap = spamreader.next()[1]
                paramExpander.cmap_cb.SetValue(cmap)
                self.model.cmap = cmap
                blur = bool_from_str(spamreader.next()[1])
                paramExpander.AdaptiveGaussianBlurBtn.SetValue(blur)
                self.model.adaptive_gaussian_blur = blur
                log = bool_from_str(spamreader.next()[1])
                paramExpander.log_CheckBox.SetValue(log)
                self.UpdateLogScale(log)
                fraction = spamreader.next()[1]
                paramExpander.fraction_cb.SetValue(fraction)
                self.UpdateFraction(fraction)
                thread = bool_from_str(spamreader.next()[1])
                paramExpander.ThreadCheckBox.SetValue(thread)
                self.use_separated_thread = thread
        except:
            None

    def BuildSource(self, event):
        WID = wid()
        self.model.update_cameras_dict(self.widgets[WID.RF_TAB].GetWindowSizesDict())
        cam = self.model.cam_dict["los"].copy()
        # extend loading by 2 AMR level:
        SelectParamDlg = SelectParametersDialog(cam.get_required_resolution() + 2)
        SelectParamDlg.ShowModal()
        if SelectParamDlg.load:
            lev = int(SelectParamDlg.levTextCtrl.GetValue())
            cam.map_max_size = 2 ** lev * max(cam.region_size) - 10
            # load every fields:
            field_list = []
            if SelectParamDlg.rhoCheckBox.GetValue():
                field_list.append("rho")
            if SelectParamDlg.PCheckBox.GetValue():
                field_list.append("P")
            if SelectParamDlg.velCheckBox.GetValue():
                field_list.append("vel")
            ngrid_max = int(float(SelectParamDlg.ngrid_maxTextCtrl.GetValue()))
            if self.use_separated_thread:
                """call BuildSource() within a dedicated thread"""

                class WorkerThread(Thread):
                    """Worker Thread Class."""

                    def __init__(self, model, cam, fl, ngrid_max):
                        super(WorkerThread, self).__init__()
                        self.model = model
                        self.cam, self.fl = cam, fl
                        self.ngrid_max = ngrid_max
                        self.start()

                    def run(self):
                        # This is the code executing in the new thread.
                        self.model.BuildSource(self.cam,
                                               self.fl, self.ngrid_max)

                # Trigger the worker thread
                WorkerThread(self.model, cam, field_list, ngrid_max)
            else:
                self.model.BuildSource(cam, field_list, ngrid_max)
        SelectParamDlg.Destroy()

    def GetCurrentMapEngine(self):
        """
        Get the good map engine corresponding to the current tab
        """
        tab_name = self.active_tab
        if tab_name is None:  # regionFinder tab case
            me_name = self.model.regionFinderMapType
        else:
            me_name = tab_name
        if me_name not in self.model.map_engine_dict:
            self.model.map_engine_dict[me_name] = map_engine(self.model, me_name, self.model.ro, self.model.source)
        return self.model.map_engine_dict[me_name]

    def UpdateLogScale(self, log_scale=None):
        """
        Update Log Scale in the appropriate MapEngine
        If called with log_scale=None, update the GUI view
        """
        tab_name = self.active_tab
        me = self.GetCurrentMapEngine()
        if log_scale is None:
            # we don't know the log_scale value
            if tab_name in self.log_scale_dict:
                # first we check into the tab dictionary
                log_scale = self.log_scale_dict[tab_name]
            else:
                # else we take the map engine default value
                log_scale = me.log_sensitive
        paramExpander = self.widgets[wid().EXPANDER_PARAMETERS]
        paramExpander.log_CheckBox.SetValue(log_scale)
        self.log_scale_dict[tab_name] = log_scale
        me.log_sensitive = log_scale

    def UpdateFraction(self, fraction=None):
        """
        Update fraction in the appropriate MapEngine
        If called with fraction=None, update the GUI view
        """
        tab_name = self.active_tab
        me = self.GetCurrentMapEngine()
        if fraction is None:
            # we don't know the fraction value
            if tab_name in self.fraction_dict:
                # first we check into the tab dictionary
                fraction = self.fraction_dict[tab_name]
            else:
                # else we take the .pymses.cfg default value
                fraction = self.default_fraction
        paramExpander = self.widgets[wid().EXPANDER_PARAMETERS]
        me.fraction = float(fraction)
        paramExpander.fraction_cb.SetValue(fraction)
        self.fraction_dict[tab_name] = fraction
        self.default_fraction = fraction

    def UpdateImages(self, wx_event=None, do_all_processing=True, verbose=True):
        if self.use_separated_thread:
            self.UpdateImagesThread(wx_event, do_all_processing, verbose)
        else:
            self.UpdateImagesNoThread(wx_event, do_all_processing, verbose)
            if self.active_tab is None:
                rf_view = self.widgets[wid().RF_TAB]
                rf_view.RefreshWholeTab()

    def UpdateImagesThread(self, wx_event, do_all_processing, verbose):
        """call UpdateImages() within a dedicated thread"""
        # Define notification event for thread completion
        EVT_RESULT_ID = wx.NewId()

        def EVT_RESULT(win, func):
            """Define Result Event."""
            win.Connect(-1, -1, EVT_RESULT_ID, func)

        class ResultEvent(wx.PyEvent):
            """Simple event to carry arbitrary result data."""

            def __init__(self, data):
                """Init Result Event."""
                super(ResultEvent, self).__init__()
                self.SetEventType(EVT_RESULT_ID)

                # self.data = data

        rf_view = self.widgets[wid().RF_TAB]
        # Set up event handler for any worker thread results
        EVT_RESULT(rf_view, rf_view.RefreshWholeTab)

        class WorkerThread(Thread):
            """Worker Thread Class."""

            def __init__(self, controller, do_all_processing, verbose):
                super(WorkerThread, self).__init__()
                self.controller = controller
                self.do_all_processing, self.verbose = do_all_processing, verbose
                self.start()

            def run(self):
                # This is the code executing in the new thread.
                self.controller.UpdateImagesNoThread(None, self.do_all_processing, self.verbose)
                if self.controller.active_tab is None:
                    rf_view = self.controller.widgets[wid().RF_TAB]
                    wx.PostEvent(rf_view, ResultEvent(None))

        # Trigger the worker thread
        WorkerThread(self, do_all_processing, verbose)

    def UpdateImagesNoThread(self, wx_event, do_all_processing, verbose):
        if self.active_tab is None:
            if "regionFinder" in self.model.image_computed or do_all_processing:
                # this test avoid a bug when user try to change the color
                # while image is not already computed
                self.RegionFinderButtonClick(do_all_processing, verbose)
        elif self.active_tab in self.model.image_computed or do_all_processing:
            # this test avoid a bug when user try to change the color
            # while image is not already computed
            self.TabMapWindowUpdateClick(do_all_processing=do_all_processing, verbose=verbose)

    def SetRegionFinderHelpMessage(self):
        self.gui.statusBar.SetStatusText(
            " Left click on two axis views to define a zoom box center (right click = cancel)", 0)
        self.gui.statusBar.SetStatusText(" Mouse wheel = zoom box size", 1)
