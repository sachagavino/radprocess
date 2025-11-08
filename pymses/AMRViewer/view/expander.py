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
import numpy
from PIL import Image, ImageDraw
from .widget_ids import WidgetIds as wid
from .utils import get_icon
from .axes3d import Axes3D
from .glass import MagGlass
from .map_tabs import GenericTabMapWindow
from pymses import rcConfig as rcConfiguration
from pymses.utils import misc
from pymses.analysis.visualization.fft_projection.convolution_kernels import ConvolKernel
from pymses.analysis.visualization.raytracing.ray_trace import ray_trace_cartesian
from pymses.analysis.visualization.image_plot_utils import *
from threading import *
import pylab as m

jet_black_cdict = {
    'red': ((0., 0., 0.), (0.30, 0.000, 0.000), (0.40, 0.2778, 0.2778), (0.52, 0.2778, 0.2778), (0.64, 1.000, 1.000),
            (0.76, 1.000, 1.000), (0.88, 0.944, 0.944), (0.98, 0.500, 0.500), (1., 1., 1.)),
    'green': ((0., 0., 0.), (0.10, 0.000, 0.000), (0.25, 0.389, 0.389), (0.32, 0.833, 0.833), (0.40, 1.000, 1.000),
              (0.52, 1.000, 1.000), (0.64, 0.803, 0.803), (0.76, 0.389, 0.389), (0.88, 0.000, 0.000), (1., 0., 0.)),
    'blue': ((0., 0.00, 0.00), (0.001, 0., 0.), (0.07, 0.500, 0.500), (0.12, 0.900, 0.900), (0.23, 1.000, 1.000),
             (0.28, 1.000, 1.000), (0.40, 0.722, 0.722), (0.52, 0.2778, 0.2778), (0.64, 0.000, 0.000), (1., 0., 0.))
}


class Expander(wx.CollapsiblePane):
    def __init__(self, parent, id, title, model, controller, is_collapsed=True):
        super(Expander, self).__init__(parent, id, label=title, style=wx.CP_DEFAULT_STYLE | wx.CP_NO_TLW_RESIZE)

        self.model = model
        self.controller = controller

        pane = self.GetPane()
        sizer = wx.BoxSizer(wx.VERTICAL)
        pane.SetSizer(sizer)

        box = wx.StaticBox(pane, wx.ID_ANY)
        boxsizer = wx.StaticBoxSizer(box, wx.VERTICAL)
        inpane = wx.Panel(pane, wx.ID_ANY)
        self.insizer = wx.BoxSizer(wx.VERTICAL)
        inpane.SetSizer(self.insizer)

        self.make_content(inpane)

        boxsizer.Add(inpane, 0, wx.ALL | wx.EXPAND, 2)
        sizer.Add(boxsizer, 0, wx.EXPAND | wx.LEFT, 10)
        self.Collapse(is_collapsed)

    def add_descr_label(self, parent, text):
        label = wx.StaticText(parent, wx.ID_ANY, "%s :" % text)
        self.insizer.Add(label, 0, wx.ALIGN_LEFT)

    def add_centered_item(self, item, is_expanded=False):
        if is_expanded:
            self.insizer.Add(item, 0, wx.LEFT | wx.EXPAND, 5)
        else:
            self.insizer.Add(item, 0, wx.LEFT | wx.ALIGN_CENTER_HORIZONTAL, 5)

    def make_content(self, parent, sizer):
        raise NotImplementedError()

    def reset(self):
        raise NotImplementedError()


class ParametersExpander(Expander):
    def __init__(self, parent, model, controller):
        WID = wid()
        super(ParametersExpander, self).__init__(parent, WID.EXPANDER_PARAMETERS, "Parameters", model, controller,
                                                 is_collapsed=True)

    def make_content(self, parent):
        WID = wid()

        self.add_descr_label(parent, "Region finder map")
        p_one = wx.Panel(parent, WID.ANY)
        self.regionFinderMapType_cb = wx.ComboBox(p_one, WID.ANY, size=(120, -1), style=wx.CB_READONLY)
        self.listregionFinderMapType = ["density", "temperature", "velocity", "Sigma", "dm_particles",
                                        "stars_particles"]
        self.regionFinderMapType_cb.SetItems(self.listregionFinderMapType)
        self.regionFinderMapType_cb.SetValue("density")
        self.model.regionFinderMapType = "density"
        self.regionFinderMapType_cb.SetToolTip(wx.ToolTip("Change the region finder map type"))
        self.add_centered_item(p_one)

        self.add_descr_label(parent, "Multiprocessing")
        p_two = wx.Panel(parent, WID.ANY)
        bs = wx.BoxSizer(wx.HORIZONTAL)
        p_two.SetSizer(bs)
        self.Multiprocessing_cb = wx.ComboBox(p_two, WID.ANY, size=(80, -1), style=wx.CB_READONLY)
        self.listMultiprocessing_cb = ["1", "2", "4", "8", "16"]
        self.Multiprocessing_cb.SetItems(self.listMultiprocessing_cb)
        self.Multiprocessing_cb.SetValue(str(rcConfiguration.multiprocessing_max_nproc))
        self.Multiprocessing_cb.SetToolTip(wx.ToolTip("Limit the number of process created with multiprocessing"))
        self.MultiprocessingBtn = wx.lib.buttons.GenBitmapToggleButton(p_two, WID.ANY,
                                                                       wx.Bitmap(get_icon("memory.png")))
        self.MultiprocessingBtn.SetToolTip(wx.ToolTip(
            "Use multiprocessing to draw the 3 region finder map in parallel, USE ONLY IF YOU HAVE ENOUGH RAM MEMORY"))
        self.MultiprocessingBtn.SetToggle(False)
        self.model.multiprocessing = self.MultiprocessingBtn.GetToggle()
        bs.Add(self.Multiprocessing_cb, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 5)
        bs.Add(self.MultiprocessingBtn, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 5)
        self.add_centered_item(p_two)

        p_resolution = wx.Panel(parent, WID.ANY)
        bs = wx.BoxSizer(wx.HORIZONTAL)
        p_resolution.SetSizer(bs)
        p_resolution_text = wx.StaticText(p_resolution, wx.ID_ANY, "Resolution : ")
        resolution_toolTip_text = "This multiply computed images resolution by this lower resolution factor to get views faster"
        p_resolution_text.SetToolTip(wx.ToolTip(resolution_toolTip_text))
        bs.Add(p_resolution_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 5)
        self.lower_resolution_cb = wx.ComboBox(p_resolution, WID.ANY, size=(60, -1), style=wx.CB_READONLY)
        self.lower_resolution_list = ["1", "0.8", "0.5", "0.3", "0.1"]
        self.lower_resolution_cb.SetItems(self.lower_resolution_list)
        self.lower_resolution_cb.SetValue("1")
        self.model.rf_lower_resolution = None
        self.model.tab_lower_resolution = None
        self.lower_resolution_cb.SetToolTip(wx.ToolTip(resolution_toolTip_text))
        bs.Add(self.lower_resolution_cb, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 5)
        # self.add_centered_item(p_resolution)
        self.insizer.Add(p_resolution, 0, wx.LEFT, 0)

        self.add_descr_label(parent, "FFT kernel size factor")
        p_three = wx.Panel(parent, WID.ANY)
        bs = wx.BoxSizer(wx.HORIZONTAL)
        p_three.SetSizer(bs)
        self.FFTkernelSize_cb = wx.ComboBox(p_three, WID.ANY, size=(80, -1), style=wx.CB_READONLY)
        self.listFFTkernelSizeStringValue = ["0.5", "1", "2", "4"]
        self.FFTkernelSize_cb.SetItems(self.listFFTkernelSizeStringValue)
        self.FFTkernelSize_cb.SetValue("1")
        self.FFTkernelSize_cb.SetToolTip(wx.ToolTip("Update FFT kernel ( = point) size factor"))
        self.model.FFTkernelSizeFactor = 1
        self.rememberSomeDataBtn = wx.lib.buttons.GenBitmapToggleButton(p_three, WID.ANY,
                                                                        wx.Bitmap(get_icon("memory.png")))
        self.rememberSomeDataBtn.SetToolTip(wx.ToolTip(
            "Remember some data in cache, USE ONLY IF YOU HAVE ENOUGH RAM MEMORY, and click here again to free all memory"))
        self.rememberSomeDataBtn.SetToggle(True)
        self.model.rememberSomeData = self.rememberSomeDataBtn.GetToggle()
        bs.Add(self.FFTkernelSize_cb, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 5)
        bs.Add(self.rememberSomeDataBtn, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 5)
        self.add_centered_item(p_three)
        self.random_shiftCheckBox = wx.CheckBox(parent, -1, 'random shift', (20, 100))
        self.random_shiftCheckBox.SetToolTip(
            wx.ToolTip("Add a random shift to point positions to avoid grid effect on FFT splatting images"))
        self.model.random_shift = False
        self.add_centered_item(self.random_shiftCheckBox)

        self.add_descr_label(parent, "Image postprocessing")
        p_four = wx.Panel(parent, WID.ANY)
        bs = wx.BoxSizer(wx.HORIZONTAL)
        p_four.SetSizer(bs)
        self.cmap_cb = wx.ComboBox(p_four, WID.ANY, size=(120, -1), style=wx.CB_READONLY)
        self.list_cmap = ["jet", "jet_black", "Greys_r", "Reds", "Greens", "Blues", "Reds_r", "Greens_r", "Blues_r",
                          "Spectral", "RdBu", "RdBu_r"]
        self.cmap_cb.SetItems(self.list_cmap)
        self.cmap_cb.SetValue("jet")
        self.cmap_cb.SetToolTip(wx.ToolTip("Change the colormap used"))
        self.AdaptiveGaussianBlurBtn = wx.lib.buttons.GenBitmapToggleButton(p_four, WID.ANY,
                                                                            wx.Bitmap(get_icon("gaussian_blur.png")))
        self.AdaptiveGaussianBlurBtn.SetToolTip(
            wx.ToolTip("Experimental : compute local image resolution and apply an adaptive gaussian blur"))
        self.AdaptiveGaussianBlurBtn.SetToggle(False)
        self.model.adaptive_gaussian_blur = self.AdaptiveGaussianBlurBtn.GetToggle()
        bs.Add(self.cmap_cb, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 5)
        bs.Add(self.AdaptiveGaussianBlurBtn, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 5)
        self.add_centered_item(p_four)
        self.log_CheckBox = wx.CheckBox(parent, -1, 'log scale', (20, 100))
        self.log_CheckBox.SetToolTip(wx.ToolTip("Turn log scale On or Off"))
        self.log_CheckBox.SetValue(True)
        self.insizer.Add(self.log_CheckBox, 0, wx.ALIGN_LEFT)
        # self.add_centered_item(self.log_CheckBox)

        p_fraction = wx.Panel(parent, WID.ANY)
        bs = wx.BoxSizer(wx.HORIZONTAL)
        p_fraction.SetSizer(bs)
        p_fraction_text = wx.StaticText(p_fraction, wx.ID_ANY, "Fraction : ")
        fraction_toolTip_text = "fraction of the total map values below the min. map range (in percent)"
        p_fraction_text.SetToolTip(wx.ToolTip(fraction_toolTip_text))
        bs.Add(p_fraction_text, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 5)
        self.fraction_cb = wx.ComboBox(p_fraction, WID.ANY, size=(80, -1), style=wx.CB_READONLY)
        self.fractionStringValue = ["0.0", "0.01", "0.1", "0.2"]
        self.fraction_cb.SetItems(self.fractionStringValue)
        self.fraction_cb.SetValue("0.1")
        self.fraction_cb.SetToolTip(wx.ToolTip(fraction_toolTip_text))
        bs.Add(self.fraction_cb, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT, 5)
        self.insizer.Add(p_fraction, 0, wx.LEFT, 0)

        # separate Thread
        self.insizer.Add(wx.StaticText(parent, wx.ID_ANY, "         -------------------------"), 0, wx.ALIGN_LEFT)
        self.controller.use_separated_thread = False
        self.ThreadCheckBox = wx.CheckBox(parent, -1, 'Separate GUI thread', (40, 100))
        self.ThreadCheckBox.SetToolTip(
            wx.ToolTip("Use a separated from GUI thread : don't use this if the GUI crash !"))
        self.insizer.Add(self.ThreadCheckBox, 0, wx.RIGHT | wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 5)

        # Bind the button
        self.ThreadCheckBox.Bind(wx.EVT_CHECKBOX, self.update_ThreadCheckBox)
        self.lower_resolution_cb.Bind(wx.EVT_COMBOBOX, self.update_lower_resolution_cb)
        self.FFTkernelSize_cb.Bind(wx.EVT_COMBOBOX, self.updateFFTkernelSize)
        self.random_shiftCheckBox.Bind(wx.EVT_CHECKBOX, self.updateRandom_shift)
        self.rememberSomeDataBtn.Bind(wx.EVT_BUTTON, self.rememberSomeDataMethod)
        self.Multiprocessing_cb.Bind(wx.EVT_COMBOBOX, self.updateMultiprocessingLimit)
        self.MultiprocessingBtn.Bind(wx.EVT_BUTTON, self.MultiprocessingBtnMethod)
        self.AdaptiveGaussianBlurBtn.Bind(wx.EVT_BUTTON, self.AdaptiveGaussianBlurBtnMethod)
        self.regionFinderMapType_cb.Bind(wx.EVT_COMBOBOX, self.updateRegionFinderMapType)
        self.cmap_cb.Bind(wx.EVT_COMBOBOX, self.update_cmap)
        self.log_CheckBox.Bind(wx.EVT_CHECKBOX, self.update_log)
        self.fraction_cb.Bind(wx.EVT_COMBOBOX, self.update_fraction)

    def update_lower_resolution_cb(self, event=None, resolution=None, save_CSV=True):
        if resolution is None:
            resolution = self.lower_resolution_list[event.GetSelection()]
        if resolution == "1":
            self.model.rf_lower_resolution = None
            self.model.tab_lower_resolution = None
        else:
            nxy_list = self.controller.widgets[wid().RF_TAB].GetWindowSizesDict()
            self.model.update_cameras_dict(nxy_list)
            self.model.rf_lower_resolution = int(float(resolution) \
                                                 * self.model.cam_dict['los'].map_max_size)
            self.model.tab_lower_resolution = int(self.model.rf_lower_resolution * 2.1)
        if event is not None:
            self.controller.UpdateImages()
        if save_CSV:
            self.controller.SaveConfigInCSV()

    def update_ThreadCheckBox(self, event):
        self.controller.use_separated_thread = self.ThreadCheckBox.GetValue()
        self.controller.SaveConfigInCSV()

    def updateRandom_shift(self, event):
        self.model.random_shift = self.random_shiftCheckBox.GetValue()
        print("random_shift =", self.model.random_shift)
        self.controller.SaveConfigInCSV()

    def updateFFTkernelSize(self, event=None, size=None, save_CSV=True):
        if size is None:
            size = self.listFFTkernelSizeStringValue[event.GetSelection()]
        self.model.FFTkernelSizeFactor = float(size)
        if event is not None:
            self.controller.UpdateImages()
        if save_CSV:
            self.controller.SaveConfigInCSV()

    def rememberSomeDataMethod(self, event):
        self.rememberSomeDataBtn.SetToggle(self.rememberSomeDataBtn.GetToggle())
        self.model.rememberSomeData = self.rememberSomeDataBtn.GetToggle()
        if not self.model.rememberSomeData:
            self.model.freeTheMemory()
        print("Remember some data =", self.model.rememberSomeData)
        self.controller.SaveConfigInCSV()

    def updateMultiprocessingLimit(self, event):
        pass
        # nproc = int(self.listMultiprocessing_cb[event.GetSelection()])
        # self.controller.SaveConfigInCSV()

    def MultiprocessingBtnMethod(self, event):
        self.MultiprocessingBtn.SetToggle(self.MultiprocessingBtn.GetToggle())
        self.model.multiprocessing = self.MultiprocessingBtn.GetToggle()
        print("Region Finder Multiprocessing =", self.model.multiprocessing)
        self.controller.SaveConfigInCSV()

    def AdaptiveGaussianBlurBtnMethod(self, event):
        self.AdaptiveGaussianBlurBtn.SetToggle(self.AdaptiveGaussianBlurBtn.GetToggle())
        self.model.adaptive_gaussian_blur = self.AdaptiveGaussianBlurBtn.GetToggle()
        self.controller.UpdateImages(do_all_processing=False, verbose=False)
        self.controller.SaveConfigInCSV()

    def updateRegionFinderMapType(self, event):
        self.model.regionFinderMapType = self.listregionFinderMapType[event.GetSelection()]
        self.model.me = None
        print("Region finder map :", self.model.regionFinderMapType)
        if self.model.regionFinderMapType == "velocity":
            self.cmap_cb.SetValue("RdBu_r")
            self.model.cmap = "RdBu_r"
            self.controller.log_scale_dict[None] = False
        else:
            self.cmap_cb.SetValue("jet")
            self.model.cmap = "jet"
            self.controller.log_scale_dict[None] = True
        self.controller.UpdateLogScale()
        self.controller.SaveConfigInCSV()

    def update_cmap(self, event):
        self.model.cmap = self.list_cmap[event.GetSelection()]
        if self.model.cmap == "jet_black":
            # generate the colormap with 1e5 interpolated values
            self.model.cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', jet_black_cdict, 1e5)
        self.controller.UpdateImages(do_all_processing=False, verbose=False)
        self.controller.SaveConfigInCSV()

    def update_log(self, event):
        log_scale = self.log_CheckBox.GetValue()
        print("log_scale =", log_scale)
        self.controller.UpdateLogScale(log_scale)
        self.controller.UpdateImages(do_all_processing=False, verbose=False)
        self.controller.SaveConfigInCSV()

    def update_fraction(self, event):
        fraction = self.fraction_cb.GetValue()
        print("fraction =", fraction)
        self.controller.UpdateFraction(fraction)
        self.controller.UpdateImages(do_all_processing=False, verbose=False)
        self.controller.SaveConfigInCSV()

    def reset(self):
        self.FFTkernelSize_cb.SetValue("1")


class LineOfSightExpander(Expander):
    def __init__(self, parent, model, controller):
        WID = wid()
        super(LineOfSightExpander, self).__init__(parent, WID.EXPANDER_LOS, "Line-of-sight axis", model, controller,
                                                  is_collapsed=False)

    def make_content(self, parent):
        WID = wid()
        self.add_descr_label(parent, "3D-view")
        self.a3d = Axes3D(parent)
        self.add_centered_item(self.a3d)
        self.a3d.Bind(wx.EVT_MOUSE_EVENTS, self.controller.A3dOnMouse)
        self.a3d.Bind(wx.EVT_KEY_DOWN, self.controller.A3dOnKeyPress)
        self.a3d.Bind(wx.EVT_KEY_UP, self.controller.A3dOnKeyRelease)

        self.add_descr_label(parent, "Axis shortcuts")
        axis_up_toolbar = wx.ToolBar(parent, WID.ANY, style=wx.TB_HORIZONTAL | wx.NO_BORDER)
        axis_up_toolbar.AddSimpleTool(WID.LOS_TOOL_X, wx.Bitmap(get_icon("x-axis.png")), 'View along the X axis')
        axis_up_toolbar.AddSeparator()
        axis_up_toolbar.AddSimpleTool(WID.LOS_TOOL_Y, wx.Bitmap(get_icon("y-axis.png")), 'View along the Y axis')
        axis_up_toolbar.AddSeparator()
        axis_up_toolbar.AddSimpleTool(WID.LOS_TOOL_Z, wx.Bitmap(get_icon("z-axis.png")), 'View along the Z axis')
        axis_up_toolbar.AddSeparator()
        axis_up_toolbar.AddSimpleTool(WID.LOS_TOOL_ROT_LEFT, wx.Bitmap(get_icon("rotate_left.png")),
                                      'Rotate the axes to the left')
        axis_up_toolbar.AddSeparator()
        axis_up_toolbar.AddSimpleTool(WID.LOS_TOOL_ROT_TOP, wx.Bitmap(get_icon("rotate_top.png")),
                                      'Rotate the axes to the top')
        axis_up_toolbar.Realize()
        axis_down_toolbar = wx.ToolBar(parent, WID.ANY, style=wx.TB_HORIZONTAL | wx.NO_BORDER)
        axis_down_toolbar.AddSimpleTool(WID.LOS_TOOL_MX, wx.Bitmap(get_icon("x-axis.png")), 'View along the -X axis')
        axis_down_toolbar.AddSeparator()
        axis_down_toolbar.AddSimpleTool(WID.LOS_TOOL_MY, wx.Bitmap(get_icon("y-axis.png")), 'View along the -Y axis')
        axis_down_toolbar.AddSeparator()
        axis_down_toolbar.AddSimpleTool(WID.LOS_TOOL_MZ, wx.Bitmap(get_icon("z-axis.png")), 'View along the -Z axis')
        axis_down_toolbar.AddSeparator()
        axis_down_toolbar.AddSimpleTool(WID.LOS_TOOL_ROT_RIGHT, wx.Bitmap(get_icon("rotate_right.png")),
                                        'Rotate the axes to the right')
        axis_down_toolbar.AddSeparator()
        axis_down_toolbar.AddSimpleTool(WID.LOS_TOOL_ROT_BOTTOM, wx.Bitmap(get_icon("rotate_bottom.png")),
                                        'Rotate the axes to the bottom')
        axis_down_toolbar.Realize()
        self.add_centered_item(axis_up_toolbar, is_expanded=True)
        self.add_centered_item(axis_down_toolbar, is_expanded=True)

        self.autorefreshCheckBox = wx.CheckBox(parent, -1, 'Auto-update on release', (20, 100))
        self.autorefreshCheckBox.SetToolTip(
            wx.ToolTip("Automatically update images when the mouse is released on 3D-view"))
        self.controller.A3d_autorefresh = False
        self.insizer.Add(self.autorefreshCheckBox, 0, wx.LEFT, 0)

        self.buildPreviewCubeCheckBox = wx.CheckBox(parent, -1, 'Build new preview cube', (20, 100))
        self.buildPreviewCubeCheckBox.SetToolTip(
            wx.ToolTip("Use amr2cube to create a low resolution preview for 3D-view"))
        self.controller.PreviewCube = False
        self.insizer.Add(self.buildPreviewCubeCheckBox, 0, wx.LEFT, 0)

        self.add_descr_label(parent, "Line-of-sight coordinates")
        coord_panel = wx.Panel(parent, WID.ANY)
        grid_sizer = wx.FlexGridSizer(3, 2)
        coord_panel.SetSizer(grid_sizer)
        los_x_label = wx.StaticText(coord_panel, WID.ANY, "x :")
        los_y_label = wx.StaticText(coord_panel, WID.ANY, "y :")
        los_z_label = wx.StaticText(coord_panel, WID.ANY, "z :")
        self.los_x = wx.TextCtrl(coord_panel, WID.ANY, style=wx.TE_READONLY, size=(80, -1))
        self.los_y = wx.TextCtrl(coord_panel, WID.ANY, style=wx.TE_READONLY, size=(80, -1))
        self.los_z = wx.TextCtrl(coord_panel, WID.ANY, style=wx.TE_READONLY, size=(80, -1))
        grid_sizer.AddMany([(los_x_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT),
                            (self.los_x, 0),
                            (los_y_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT),
                            (self.los_y, 0),
                            (los_z_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT),
                            (self.los_z, 0)])
        self.add_centered_item(coord_panel)

        # Binding toolbar events to controller method
        axis_up_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_X)
        axis_down_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_MX)
        axis_up_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_Y)
        axis_down_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_MY)
        axis_up_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_Z)
        axis_down_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_MZ)
        axis_up_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_ROT_LEFT)
        axis_down_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_ROT_RIGHT)
        axis_up_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_ROT_TOP)
        axis_down_toolbar.Bind(wx.EVT_TOOL, self.controller.OnShortcutAxis, id=WID.LOS_TOOL_ROT_BOTTOM)
        self.autorefreshCheckBox.Bind(wx.EVT_CHECKBOX, self.updateAutorefresh)
        self.buildPreviewCubeCheckBox.Bind(wx.EVT_CHECKBOX, self.updatePreviewCube)

    def RefreshAxes3D(self):
        """
        Line-of-sight expander content refresh method
        """
        # Axes3D refresh
        self.a3d.SetUVPoints(self.model.GetUVPoints())
        self.a3d.Refresh()

        # Line-of-sight x/y/z coordinates update
        los = self.model.GetLosAxis()
        self.los_x.ChangeValue("%.6f" % los[0])
        self.los_y.ChangeValue("%.6f" % los[1])
        self.los_z.ChangeValue("%.6f" % los[2])

        # real time preview
        # from time import time
        if self.controller.PreviewCube:
            # t0 = time()
            mms = numpy.min(self.model.preview_cube.shape)
            cam = self.model.get_los_camera(mms)
            (ray_vectors, ray_origins, ray_lengths) = cam.get_rays()
            map = ray_trace_cartesian(self.model.preview_cube, ray_origins, ray_vectors, ray_lengths,
                                      cam.get_bounding_box(), mms)
            map = map.reshape(mms, mms)
            map = map.transpose()

            map = apply_log_scale(map, verbose=False)
            # import pylab
            # imgPylab=pylab.imshow(map)
            # pylab.show()
            # vmin, vmax = numpy.min(map), numpy.max(map)
            # print numpy.min(map), numpy.max(map)
            # colormap = pylab.cm.get_cmap("jet")
            vmin, vmax = get_map_range(map, log_sensitive=True)
            colormap, vmin, vmax = get_Colormap("jet", False, (vmin, vmax))
            map = (map - vmin) / (vmax - vmin)
            map = numpy.asarray(colormap(map) * 255, dtype='i')
            map = map.reshape(mms * mms, 4)
            R_band = Image.new("L", (mms, mms))
            R_band.putdata(map[:, 0])
            G_band = Image.new("L", (mms, mms))
            G_band.putdata(map[:, 1])
            B_band = Image.new("L", (mms, mms))
            B_band.putdata(map[:, 2])
            im = Image.merge("RGB", (R_band, G_band, B_band)).transpose(Image.FLIP_TOP_BOTTOM)
            im = im.resize((200, 200))
            # im.show()
            image = wx.EmptyImage(im.size[0], im.size[1])
            image.SetData(im.convert("RGB").tostring())
            bitmap = wx.BitmapFromImage(image)
            self.a3d.SetPreviewImage(bitmap)
            # print "compute image :", time()-t0
            self.a3d.OnPaint(None)  # make sure the new image is painted

    def GetAxes3d(self):
        return self.a3d

    def updateAutorefresh(self, event):
        self.controller.A3d_autorefresh = self.autorefreshCheckBox.GetValue()

    def updatePreviewCube(self, event):
        self.controller.PreviewCube = self.buildPreviewCubeCheckBox.GetValue()
        if self.controller.PreviewCube:
            if self.controller.use_separated_thread:
                # Define notification event for thread completion
                EVT_RESULT_ID = wx.NewId()

                def EVT_RESULT(win, func):
                    win.Connect(-1, -1, EVT_RESULT_ID, func)

                class ResultEvent(wx.PyEvent):
                    def __init__(self):
                        super(ResultEvent, self).__init__()
                        self.SetEventType(EVT_RESULT_ID)

                # Set up event handler for any worker thread results
                EVT_RESULT(self, self.RefreshAxes3D)

                """call BuildPreviewCube() within a dedicated thread"""

                class WorkerThread(Thread):
                    def __init__(self, parent_thread):
                        super(WorkerThread, self).__init__()
                        self.parent_thread = parent_thread

                    def run(self):
                        self.parent_thread.model.BuildPreviewCube()
                        wx.PostEvent(self.parent_thread, ResultEvent())

                WorkerThread(self).start()
            else:
                self.model.BuildPreviewCube()
        else:
            self.a3d.SetPreviewImage(None)
        self.RefreshAxes3D()

    def reset(self):
        self.a3d.SetCtrlDown(False)
        self.a3d.SetClickDown(False)
        self.a3d.SetMouseIn(False)
        self.a3d.SetCursorPosition(None, None)
        self.a3d.UpdateCursor()
        self.RefreshAxes3D()


class MagGlassExpander(Expander):
    def __init__(self, parent, model, controller):
        WID = wid()
        super(MagGlassExpander, self).__init__(parent, WID.EXPANDER_GLASS, "Magnifier", model, controller,
                                               is_collapsed=True)

    def make_content(self, parent):
        WID = wid()
        self.add_descr_label(parent, "Magnifying glass")
        self.glass = MagGlass(parent)
        self.add_centered_item(self.glass)

        self.add_descr_label(parent, "Value")
        self.glass_value = wx.StaticText(parent, WID.GLASS_VAL, "        ---        ")
        self.add_centered_item(self.glass_value)

        self.add_descr_label(parent, "Magnification factor")
        self.mag_factor_list = [2, 4, 5, 8, 10, 20]
        self.mag_factor_cb = wx.ComboBox(parent, WID.ANY, size=(80, -1), style=wx.CB_READONLY)
        self.mag_factor_cb.AppendItems(["x %i" % f for f in self.mag_factor_list])
        self.mag_factor_cb.Select(3)
        self.mag_factor_cb.Bind(wx.EVT_COMBOBOX, self.on_change_mag_factor)
        self.add_centered_item(self.mag_factor_cb)

    def reset(self):
        pass

    def SetGlassImage(self, pil_image):
        if pil_image is None:
            self.glass.bitmap = None
        else:
            nx, ny = self.glass.GetSize()
            img = wx.EmptyImage(nx, ny)
            rpil = pil_image.resize((nx, ny))
            img.SetData(rpil.convert("RGB").tostring())
            img.SetAlphaData(rpil.tostring()[3::4])
            self.glass.bitmap = img.ConvertToBitmap()
        self.glass.Refresh()

    def RefreshPhysicalValue(self, label):
        """
        """
        self.glass_value.SetLabel(label)

    def on_change_mag_factor(self, event):
        self.set_mag_factor(event.GetSelection())

    def set_mag_factor(self, index):
        mgf = self.mag_factor_list[index]
        nx, ny = self.glass.GetSize()
        GenericTabMapWindow.mag_box_size = nx / mgf
        MagGlass.mag_factor = mgf


__all__ = ["ParametersExpander", "LineOfSightExpander", "MagGlassExpander"]
