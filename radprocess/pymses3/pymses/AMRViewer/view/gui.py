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
from .about import open_about_box
from .left_panel import LeftPanel
from .notebook import Notebook
from .dialogs import *
from ..ctrl.controller import Controller
from .utils import get_icon
import os


class GUI(wx.App):
    def OnInit(self):
        frame = MainWindow(None, title='AMRViewer')
        self.SetTopWindow(frame)
        return True

    def run(self):
        self.MainLoop()


class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        """
        """
        WID = wid()
        super(MainWindow, self).__init__(parent, WID.ANY, title=title)

        # Controller
        self.controller = Controller(self)

        # Model
        self.model = self.controller.model

        # Setting up the menubar
        filemenu = wx.Menu()
        menuOpen = filemenu.Append(WID.MENU_OPEN, "&Open RAMSES simulation", " Open a RAMSES outputs directory")
        menuExit = filemenu.Append(WID.MENU_QUIT, "&Quit", " Terminate the AMRViewer")

        cameramenu = wx.Menu()
        menuResetCam = cameramenu.Append(WID.TOOL_RESET, "&Reset view", " Reset view")
        menuLoadCam = cameramenu.Append(WID.MENU_LOADCAM, "&Load camera...", " Load a camera from a file")
        menuSaveCam = cameramenu.Append(WID.MENU_SAVECAM, "&Save camera into file...", " Save the camera into a file")

        self.tabdisplaymenu = wx.Menu()
        for tab_index in range(len(WID.TAB_NAME_LIST)):
            tab_name = WID.TAB_NAME_LIST[tab_index]
            descr = WID.TAB_DESCR_LIST[tab_name]
            self.tabdisplaymenu.Append(WID.TAB_MENU_IDS[tab_name], "%s" % descr, " Show/hide tab : '%s'" % descr,
                                       wx.ITEM_CHECK)
            if tab_name != "hdf5" and tab_name != "custom":
                self.tabdisplaymenu.Check(WID.TAB_MENU_IDS[tab_name], True)

        helpmenu = wx.Menu()
        menuAbout = helpmenu.Append(WID.MENU_ABOUT, "&About", " Information on AMRViewer")

        menuBar = wx.MenuBar()
        menuBar.Append(filemenu, "&File")  # Adding the "filemenu" to the MenuBar
        menuBar.Append(cameramenu, "&Camera")  # Adding the "cameramenu" to the MenuBar
        menuBar.Append(self.tabdisplaymenu, "&Display")  # Adding the "tabdisplaymenu" to the MenuBar
        menuBar.Append(helpmenu, "&Help")  # Adding the "helpmenu" to the MenuBar
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

        # Toolbar creation + events
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.toolbar_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.toolbar_panel = wx.Panel(self, -1)
        self.toolbar_panel.SetSizer(self.toolbar_sizer)
        self.SetSizer(sizer)
        toolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL)
        toolbar.AddSimpleTool(WID.TOOL_OPEN, wx.Bitmap(get_icon("open.png")), '', 'Open a RAMSES outputs directory')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(WID.TOOL_QUIT, wx.Bitmap(get_icon("quit.png")), '', 'Terminate the AMRViewer')
        toolbar.AddSimpleTool(WID.TOOL_CLOSE, wx.Bitmap(get_icon("stop.png")), '', 'Close current tab')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(WID.TOOL_RESET, wx.Bitmap(get_icon("reset.png")), '', 'Reset view')
        toolbar.AddSimpleTool(WID.TOOL_PREVIOUS, wx.Bitmap(get_icon("previous.png")), '', 'Previous view')
        toolbar.AddSimpleTool(WID.MENU_LOADCAM, wx.Bitmap(get_icon("camera_open.png")), '', 'Load a camera from a file')
        toolbar.AddSimpleTool(WID.MENU_SAVECAM, wx.Bitmap(get_icon("camera_save.png")), '',
                              'Save the camera into a file')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(WID.TOOL_GLNEMO2, wx.Bitmap(get_icon("glnemo2.png")), '',
                              'Start glnemo2 on this view box')
        toolbar.AddSeparator()
        toolbar.AddSimpleTool(WID.TOOL_SAVE_IMAGE, wx.Bitmap(get_icon("save.png")), '', 'Save image')
        toolbar.AddSimpleTool(WID.TOOL_BUILD_SOURCE, wx.Bitmap(get_icon("load.png")), '',
                              'Build a local octree for this area')
        toolbar.AddSimpleTool(WID.TOOL_UPDATE, wx.Bitmap(get_icon("update.png")), '', 'Update map')

        toolbar.Realize()
        self.toolbar_sizer.Add(toolbar, wx.ALIGN_LEFT)

        self.hdf5_cb = wx.ComboBox(self.toolbar_panel, WID.ANY, size=(380, -1), style=wx.CB_READONLY)
        self.hdf5_cb.SetToolTip(wx.ToolTip("This box is only used to select files in hdf5 tab!"))
        self.toolbar_sizer.Add(self.hdf5_cb, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT, 5)
        self.hdf5_cb.Bind(wx.EVT_COMBOBOX, self.update_hdf5)
        self.hdf5_cb.Hide()

        self.stars_age_instensity_dimmingCheckBox = wx.CheckBox(self.toolbar_panel, -1,
                                                                'Stars intensity dimming with age', (20, 100))
        self.stars_age_instensity_dimmingCheckBox.SetToolTip(
            wx.ToolTip("Stars intensity dimming for stars more than 10 Myr old"))
        self.model.stars_age_instensity_dimming = True
        self.stars_age_instensity_dimmingCheckBox.SetValue(self.model.stars_age_instensity_dimming)
        self.toolbar_sizer.Add(self.stars_age_instensity_dimmingCheckBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT,
                               5)
        self.stars_age_instensity_dimmingCheckBox.Bind(wx.EVT_CHECKBOX, self.update_stars_age_instensity_dimming)
        self.stars_age_instensity_dimmingCheckBox.Hide()

        sizer.Add(self.toolbar_panel)

        # StatusBar
        self.statusBar = self.CreateStatusBar(2)  # A Statusbar in the bottom of the window
        self.SetStatusWidths([-5, -2])
        self.statusBar.SetStatusText(" Welcome in the amrviewer gui !")

        # Frame Content
        mainpanel = wx.Panel(self, WID.ANY, style=wx.RAISED_BORDER)
        sizer.Add(mainpanel, 1, wx.EXPAND)
        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        mainpanel.SetSizer(sizer2)
        self.lp = LeftPanel(mainpanel, self.model, self.controller)
        self.rp = Notebook(mainpanel, self.model, self.controller)
        sizer2.Add(self.lp, 0, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER, 2)
        sizer2.Add(self.rp, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER, 2)

        # Widgets
        self.widgets = {
            WID.EXPANDER_PARAMETERS: self.lp.cpp,
            WID.EXPANDER_LOS: self.lp.cpa,
            WID.EXPANDER_GLASS: self.lp.cpg,
            WID.RF_TAB: self.rp.rf_tab,
        }
        for tab_name, tab_id in WID.TAB_IDS.items():
            self.widgets[tab_id] = self.rp.map_tabs[tab_name]

        # Controller
        self.controller.SetWidgetsDict(self.widgets)

        # Controller reset
        self.controller.OnReset(None)

        # Events binding
        self.bind_events()

        self.SetAutoLayout(1)
        self.Fit()
        self.SetMinSize(self.GetSize())
        self.Center()
        self.Show()

        # show tab by default
        for tab_name, tab_menu_id in WID.TAB_MENU_IDS.items():
            if tab_name != "hdf5" and tab_name != "custom":
                tab = self.widgets[WID.TAB_IDS[tab_name]]
                tab.Show()

        # select current directory by default
        self.controller.SelectRamsesDir(abspath("./"))

    def reset(self):
        self.lp.reset()
        self.rp.reset()

    def bind_events(self):
        WID = wid()
        self.Bind(wx.EVT_MENU, self.OnOpen, id=WID.MENU_OPEN)
        self.Bind(wx.EVT_TOOL, self.OnExit, id=WID.MENU_QUIT)

        self.Bind(wx.EVT_MENU, self.OnLoadCam, id=WID.MENU_LOADCAM)
        self.Bind(wx.EVT_MENU, self.OnSaveCam, id=WID.MENU_SAVECAM)
        self.Bind(wx.EVT_MENU, self.OnGlnemo2, id=WID.TOOL_GLNEMO2)

        for tid in list(WID.TAB_MENU_IDS.values()):
            self.Bind(wx.EVT_MENU, self.OnDisplayHideTab, id=tid)

        self.Bind(wx.EVT_MENU, self.OnAbout, id=WID.MENU_ABOUT)

        self.Bind(wx.EVT_TOOL, self.OnOpen, id=WID.TOOL_OPEN)
        self.Bind(wx.EVT_MENU, self.OnExit, id=WID.TOOL_QUIT)
        self.Bind(wx.EVT_TOOL, self.controller.OnReset, id=WID.TOOL_RESET)
        self.Bind(wx.EVT_TOOL, self.controller.OnPrevious, id=WID.TOOL_PREVIOUS)
        self.Bind(wx.EVT_TOOL, self.controller.UpdateImages, id=WID.TOOL_UPDATE)
        self.Bind(wx.EVT_TOOL, self.onCloseTab, id=WID.TOOL_CLOSE)
        self.Bind(wx.EVT_TOOL, self.controller.BuildSource, id=WID.TOOL_BUILD_SOURCE)
        self.Bind(wx.EVT_TOOL, self.OnSaveImg, id=WID.TOOL_SAVE_IMAGE)

        self.Bind(wx.EVT_COLLAPSIBLEPANE_CHANGED, self.UpdateLayout)

    def UpdateLayout(self, event):
        self.Layout()

    def update_stars_age_instensity_dimming(self, event):
        self.model.stars_age_instensity_dimming = self.stars_age_instensity_dimmingCheckBox.GetValue()
        self.controller.UpdateImages()

    def OnAbout(self, event):
        open_about_box(self)

    def OnExit(self, event):
        dial = wx.MessageDialog(self, 'Are you sure you want to quit ?', 'Exit', \
                                wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION)
        if dial.ShowModal() == wx.ID_YES:
            self.Close(True)

    def OnSaveImg(self, event):
        """
        Open a CameraPNGFileDialog and let the user choose a PNG file
        in which a image is saved
        """
        if self.controller.active_tab is None and self.rp.rf_tab.windows[wid().RF_TAB_LOS_WIN].image is not None:
            selectViewDialog = SelectViewDialog()
            selectViewDialog.ShowModal()
            v = selectViewDialog.selected_view
            selectViewDialog.Destroy()

            win = CameraPNGFileDialog(self, self.controller, self.rp.rf_tab.windows[v].image)
            win.run()

        elif self.controller.active_tab is not None and \
                        self.rp.map_tabs[self.controller.active_tab].window.image is not None:

            typeDlg = FileTypeSelectionDialog(self, self.controller)
            typeDlg.ShowModal()
            is_png = typeDlg.GetStringSelection() == "PNG"
            typeDlg.Destroy()

            if is_png:
                win = CameraPNGFileDialog(self, self.controller,
                                          self.rp.map_tabs[self.controller.active_tab].window.image)
                win.run()
            else:
                fwin = FitsFileDialog(self, self.controller, self.model.GetMap(self.controller.active_tab),
                                      self.model.GetMapLengthUnit(), self.model.mapUnit[self.controller.active_tab])
                fwin.run()
        else:
            print("There is no image to save!")

    def OnOpen(self, event):
        """
        Open a RamsesSimOuputDialog
        """
        WID = wid()
        win = RamsesSimOutputDialog(self, self.controller)
        win.run()

    def OnLoadCam(self, event):
        """
        Open a CameraFileDialog and let the user choose a HDF5 file
        from which a Camera is loaded
        """
        WID = wid()
        win = CameraFileDialog(self, self.controller, True)
        win.run()

    def OnSaveCam(self, event):
        """
        Open a CameraFileDialog and let the user choose a HDF5 file
        in which a Camera is saved
        """
        WID = wid()
        win = CameraFileDialog(self, self.controller, False)
        win.run()

    def OnGlnemo2(self, event):
        """
        Start glnemo2 on this view box
        """
        WID = wid()
        from pymses.sources.ramses import filename_utils

        iout = self.model.iout_list[self.model.selected_iout_index]
        out_dir = filename_utils.output_path(self.model.ramses_dir, iout)
        d = self.rp.rf_tab.GetWindowSizesDict()
        self.model.update_cameras_dict(d)
        rs2 = self.model.cam_dict["los"].region_size[0] / 2
        xmin = self.model.cam_dict["los"].center[0] - rs2
        xmax = self.model.cam_dict["los"].center[0] + rs2
        ymin = self.model.cam_dict["los"].center[1] - rs2
        ymax = self.model.cam_dict["los"].center[1] + rs2
        zmin = self.model.cam_dict["los"].center[2] - rs2
        zmax = self.model.cam_dict["los"].center[2] + rs2
        lmax = self.model.cam_dict["los"].get_required_resolution()
        lmax = min(self.model.ro.info["levelmax"], lmax + 2)
        if self.controller.active_tab == "dm_particles":
            select = "halo"
        elif self.controller.active_tab == "stars_particles":
            select = "stars"
        else:
            select = "gas"
        glnemo2_cmd = "glnemo2 in=" + out_dir + " select=" + select + \
                      " xmin=" + str(xmin) + " xmax=" + str(xmax) + " ymin=" + \
                      str(ymin) + " ymax=" + str(ymax) + " zmin=" + str(zmin) \
                      + " zmax=" + str(zmax)
        # if self.controller.active_tab != "particles":
        glnemo2_cmd += " #lmax=" + str(lmax)
        dlg = wx.TextEntryDialog(self, 'glnemo2 command has to be available in your path !', \
                                 'Start Glnemo2 ?', glnemo2_cmd)
        if dlg.ShowModal() == wx.ID_OK:
            print("Starting command :", dlg.GetValue())
            os.system(dlg.GetValue())
        dlg.Destroy()

    # os.system(glnemo2_cmd)
    # import subprocess
    # subprocess.call(glnemo2_cmd, shell=True)

    def OnDisplayHideTab(self, event):
        WID = wid()
        tid = event.GetId()
        for tab_name, tab_menu_id in WID.TAB_MENU_IDS.items():
            if (tab_menu_id == tid):
                tab = self.widgets[WID.TAB_IDS[tab_name]]
                if tab.IsShown():
                    tab.Hide()
                else:
                    tab.Show()
                    self.controller.ChangeActiveTab(tab_name)
                    self.rp.ChangeSelection(WID.Tab_Number_for_selection[tab_name])
                    if tab_name == "hdf5":
                        win = HDF5FileDialog(self.controller)
                        win.run()
                        tab = self.controller.widgets[WID.TAB_IDS[tab_name]]
                        if tab_name in self.model.image_computed and \
                                        self.model.image_computed[tab_name] == True:
                            # A hdf5 file has been selected
                            tab.UpdateTabImage()
                            tab.RefreshWholeTab(make_bitmap=True)
                break

    def onCloseTab(self, event):
        tab_name = self.controller.active_tab
        if tab_name is None:
            print("Region finder cannot be closed !")
        else:
            tab = self.widgets[wid().TAB_IDS[tab_name]]
            tab.Hide()
            self.tabdisplaymenu.Check(wid().TAB_MENU_IDS[tab_name], False)

    def update_hdf5(self, event):
        if self.controller.active_tab == "hdf5":
            self.model.LoadmapFromHDF5(self.hdf5_file_dirname + "/" + self.list_hdf5_file[event.GetSelection()])
            WID = wid()
            tab = self.widgets[WID.TAB_IDS["hdf5"]]
            tab.UpdateTabImage()
            tab.RefreshWholeTab(make_bitmap=True)


def run():
    gui = GUI()
    gui.run()
