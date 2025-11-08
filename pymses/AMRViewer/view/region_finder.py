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
from .utils import get_icon
from .widget_ids import WidgetIds as wid
from numpy import abs, min, max, round
from .map_tabs import MapWindow
from .dialogs import EditCameraDialog
from pymses.analysis.visualization import Camera


class RegionFinderTab(wx.Panel):
    axis_name = ["los", "u", "v"]

    label_string = {"los": "line-of-sight view", "u": "u view", "v": "-v view"}

    def __init__(self, parent, model, controller):
        WID = wid()
        self.controller = controller
        self.axis_ids = {"los": WID.RF_TAB_LOS_WIN, "u": WID.RF_TAB_U_WIN, "v": WID.RF_TAB_V_WIN}
        super(RegionFinderTab, self).__init__(parent, WID.RF_TAB)
        self.model = model

        page_sizer = wx.GridSizer(2, 2)
        self.SetSizer(page_sizer)
        self.windows = {}
        for axis in RegionFinderTab.axis_name:
            id = self.axis_ids[axis]
            lab = RegionFinderTab.label_string[axis]
            b = wx.StaticBox(self, WID.ANY, lab)
            bs = wx.StaticBoxSizer(b, wx.VERTICAL)
            self.windows[id] = AxisWindow(self, id, axis)
            self.windows[id].Bind(wx.EVT_MOUSE_EVENTS, controller.RegionFinderOnMouse)
            bs.Add(self.windows[id], 1, wx.EXPAND | wx.SHAPED | wx.ALIGN_CENTER | wx.ALL, 1)

            page_sizer.Add(bs, 1, wx.EXPAND | wx.SHAPED | wx.ALIGN_CENTER | wx.ALL, 4)

        # Region params StaticBox
        box = wx.StaticBox(self, WID.ANY, "Region of interest")
        boxsizer = wx.StaticBoxSizer(box, wx.VERTICAL)

        pane = wx.Panel(self, WID.ANY)
        psizer = wx.BoxSizer(wx.VERTICAL)
        pane.SetSizer(psizer)

        # output
        self.out_dir_ctrl = wx.TextCtrl(pane, WID.ANY, style=wx.TE_READONLY, size=(250, -1))
        psizer.Add(self.out_dir_ctrl, 0, wx.EXPAND)

        output_pane = wx.Panel(pane, WID.ANY)
        bs = wx.BoxSizer(wx.HORIZONTAL)
        output_pane.SetSizer(bs)
        self.iout_cb = wx.ComboBox(output_pane, WID.ANY, size=(80, -1), style=wx.CB_READONLY)
        self.iout_cb.SetEditable(True)
        iout_cb_refbtn = wx.BitmapButton(output_pane, WID.ANY, wx.Bitmap(get_icon("update.png")))
        iout_cb_refbtn.SetToolTip(wx.ToolTip("Update output list"))
        output_text = wx.StaticText(output_pane, WID.ANY, "Output")

        bs.Add(output_text, flag=wx.ALIGN_CENTER_VERTICAL)
        bs.Add(self.iout_cb, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)
        bs.Add(iout_cb_refbtn, flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL, border=5)

        psizer.Add(output_pane, 0, wx.ALIGN_CENTER)

        # Bind button to controller method
        self.iout_cb.Bind(wx.EVT_COMBOBOX, controller.OnSelectIout)
        wx.EVT_TEXT_ENTER(output_pane, -1, controller.OnSelectIout)
        self.iout_cb.Bind(wx.EVT_KILL_FOCUS, controller.OnSelectIout)
        iout_cb_refbtn.Bind(wx.EVT_BUTTON, controller.RefreshRamsesIoutList)

        # Region center coordinates
        rclabel = wx.StaticText(pane, WID.ANY, "Target center coordinates :")
        psizer.Add(rclabel, 0, wx.ALIGN_LEFT)

        coord_panel = wx.Panel(pane, WID.ANY)
        grid_sizer = wx.FlexGridSizer(3, 4)
        coord_panel.SetSizer(grid_sizer)
        rc_x_label = wx.StaticText(coord_panel, WID.ANY, "x :")
        rc_y_label = wx.StaticText(coord_panel, WID.ANY, "y :")
        rc_z_label = wx.StaticText(coord_panel, WID.ANY, "z :")
        self.rc_x = wx.TextCtrl(coord_panel, WID.ANY, style=wx.TE_READONLY, size=(75, -1))
        self.rc_y = wx.TextCtrl(coord_panel, WID.ANY, style=wx.TE_READONLY, size=(75, -1))
        self.rc_z = wx.TextCtrl(coord_panel, WID.ANY, style=wx.TE_READONLY, size=(75, -1))
        empty_label1 = wx.StaticText(coord_panel, WID.ANY, "  ")
        empty_label2 = wx.StaticText(coord_panel, WID.ANY, "  ")
        empty_label3 = wx.StaticText(coord_panel, WID.ANY, "  ")
        empty_label4 = wx.StaticText(coord_panel, WID.ANY, "  ")
        self.edit_cam_btn = wx.Button(coord_panel, -1, "Edit", size=(50, -1))
        grid_sizer.AddMany([(rc_x_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT),
                            (self.rc_x, 0), (empty_label1, 0), (empty_label2, 0),
                            (rc_y_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT),
                            (self.rc_y, 0), (empty_label3, 0), (self.edit_cam_btn, 0),
                            (rc_z_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT),
                            (self.rc_z, 0), (empty_label4, 0)])

        self.edit_cam_btn.Bind(wx.EVT_BUTTON, self.edit_cam)

        psizer.Add(coord_panel, 0, wx.LEFT | wx.ALIGN_CENTER_HORIZONTAL, 5)

        # Region size
        rslabel = wx.StaticText(pane, WID.ANY, "Region/Zoom size :")
        psizer.Add(rslabel, 0, wx.ALIGN_LEFT)

        rsize_panel = wx.Panel(pane, WID.ANY)
        rs_sizer = wx.BoxSizer(wx.HORIZONTAL)
        rsize_panel.SetSizer(rs_sizer)
        rs_val_label = wx.StaticText(rsize_panel, WID.ANY, "  L =")
        self.rs = wx.TextCtrl(rsize_panel, WID.ANY, style=wx.TE_READONLY, size=(60, -1))
        self.rs_unit_label = wx.StaticText(rsize_panel, WID.ANY, "", size=(40, -1))
        rs_sizer.Add(rs_val_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT)
        rs_sizer.Add(self.rs, 0)
        rs_sizer.Add(self.rs_unit_label, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_LEFT)

        psizer.Add(rsize_panel, 0, wx.LEFT | wx.ALIGN_CENTER_HORIZONTAL, 5)

        # Output time
        out_time_label = wx.StaticText(pane, WID.ANY, "Output time :")
        self.out_time_val_label = wx.StaticText(pane, WID.ANY, "---")
        psizer.Add(out_time_label, 0, wx.ALIGN_LEFT)
        psizer.Add(self.out_time_val_label, 0, wx.ALIGN_CENTER)

        boxsizer.Add(pane, 0, wx.ALL | wx.EXPAND, 2)
        page_sizer.Add(boxsizer, 0, wx.ALIGN_CENTER)

    def edit_cam(self, evt):
        self.model.update_cameras_dict(self.GetWindowSizesDict())
        los_cam = self.model.cam_dict["los"]
        EditCameraDlg = EditCameraDialog(los_cam)
        EditCameraDlg.ShowModal()
        if EditCameraDlg.ok:
            center = eval(EditCameraDlg.center.GetValue())
            line_of_sight_axis = eval(EditCameraDlg.line_of_sight_axis.GetValue())
            rs = float(EditCameraDlg.region_size.GetValue())
            new_cam = Camera(center=center, line_of_sight_axis=line_of_sight_axis, region_size=[rs, rs],
                             size_unit=los_cam.size_unit, distance=(rs / 2.), far_cut_depth=(rs / 2.),
                             map_max_size=self.model.cam_dict["los"].map_max_size)
            self.model.LoadLosCamera(new_cam)
            self.controller.gui.reset()
        EditCameraDlg.Destroy()

    def reset(self):
        self.reset_tab()
        self.reset_windows()
        self.RefreshWholeTab()
        self.RefreshOutputInfo()

    def reset_tab(self):
        self.FreezeCursorPosition(u=False, v=False, w=False)
        self.MoveCursor(u=None, v=None, w=None)
        self.SetSizeSet(False)
        self.ResetZoomBox()

    def reset_windows(self):
        for win in list(self.windows.values()):
            win.reset()

    def RefreshWholeTab(self, make_bitmap=True):
        self.RefreshCrossCenterCoords()
        self.RefreshAxisWindows(make_bitmap)
        self.RefreshRegionSize()
        self.controller.SetRegionFinderHelpMessage()

    def RefreshOutputInfo(self, select_last_iout=False):
        if self.model.ramses_dir is not None:
            # directory
            self.out_dir_ctrl.ChangeValue(self.model.ramses_dir)
            # iout : display the output list in the ComboBox and select the first item
            ilist = self.model.GetvalidOutputsList()
            ilist_string = ["%i" % i for i in ilist]
            self.iout_cb.SetItems(ilist_string)
            if select_last_iout:
                self.controller.SelectRamsesOutput(len(ilist) - 1)

    def SetRamsesOutput(self):
        if self.model.ramses_dir is not None:
            self.iout_cb.Select(self.model.selected_iout_index)
        self.RefreshRegionSize()
        self.RefreshOutputTime()

    def RefreshOutputTime(self):
        try:  # self.model.ramses_dir is not None:
            time, unit = self.model.GetOutputTime()
            self.out_time_val_label.SetLabel("%f %s" % (time, unit))
        except:
            self.out_time_val_label.SetLabel("---")

    def RefreshCrossCenterCoords(self):
        cc = self.model.GetCrossCenter()
        self.rc_x.ChangeValue("%.6f" % cc[0])
        self.rc_y.ChangeValue("%.6f" % cc[1])
        self.rc_z.ChangeValue("%.6f" % cc[2])

    def RefreshRegionSize(self):
        rs, rs_unit = self.model.GetRegionSize()
        self.rs.ChangeValue("%.2f" % rs)
        self.rs_unit_label.SetLabel(" %s" % rs_unit)

    # Window methods
    def GetWindowSizesDict(self):
        d = {}
        for axis in RegionFinderTab.axis_name:
            id = self.axis_ids[axis]
            win = self.windows[id]
            d[axis] = win.nx
        return d

    def UpdateRegionFinderImages(self):
        d = self.model.GetAxisImages()
        for win in list(self.windows.values()):
            axis = win.axis
            if axis in d:
                win.SetBackgroundImage(d[axis])

    def RefreshAxisWindows(self, make_bitmap):
        for win in list(self.windows.values()):
            win.RefreshWindow(make_bitmap)

    def RefreshWindowCursors(self):
        """
        Updates the cursor state on each window
        """
        for win in list(self.windows.values()):
            win.RefreshCursor()

    @classmethod
    def IsZoomSizeIsSet(self):
        return AxisWindow.is_size_set

    @classmethod
    def SetSizeSet(self, s):
        AxisWindow.is_size_set = s

    @classmethod
    def ChangeZoomSize(self, wh, ds=0.02):
        """
        Modify the zoom size according to the wheel rotation value
        """
        dzs = ds
        if AxisWindow.L is not None:
            dzs *= AxisWindow.L
        if wh > 0:
            zs = AxisWindow.zoom_size - dzs
        else:
            zs = AxisWindow.zoom_size + dzs
        return self.check_valid_zoom_size(zs, border=False)

    @classmethod
    def InitZoomSize(self):
        """
        Try to initialize the zoom box to half the size of the window
        """
        self.check_valid_zoom_size(0.5)

    @classmethod
    def ResetZoomBox(self):
        AxisWindow.zoom_size = None
        AxisWindow.L = None

    @classmethod
    def check_valid_zoom_size(self, zs, border=True):
        # Min. value of the zoom size = 1% of the window
        if zs < 0.01:
            return False

        # zoom size required to keep the zoom box inside the window
        zsb = min([zs, 1.0 - 2.0 * max([abs(AxisWindow.cursor_u), abs(AxisWindow.cursor_v), abs(AxisWindow.cursor_w)])])

        if border:  # the zoom box must remain inside the window
            AxisWindow.zoom_size = zsb
            return False
        else:
            if zs <= zsb:  # No problem
                AxisWindow.zoom_size = zs
                if AxisWindow.L is not None:
                    AxisWindow.L = None
                    return True
                return False
            else:  # Max. value of the zoom size
                # Do something, we want to make the zoom box leave the window !
                AxisWindow.zoom_size = zs
                AxisWindow.L = 1.0 + zs - zsb
                return True

    @classmethod
    def GetZoomSize(self):
        """ Gets the current zoom size
        """
        zs = AxisWindow.zoom_size
        if AxisWindow.L is not None:
            zs = (zs + 1.0 - AxisWindow.L) * AxisWindow.L
        return zs

    @classmethod
    def MoveCursor(self, u=-1, v=-1, w=-1):
        """
        Move cursor position (u, v, w) if not set
        """
        if u != -1:
            if not AxisWindow.is_set_u:
                AxisWindow.cursor_u = u
        if v != -1:
            if not AxisWindow.is_set_v:
                AxisWindow.cursor_v = v
        if w != -1:
            if not AxisWindow.is_set_w:
                AxisWindow.cursor_w = w

    @classmethod
    def GetCursorUVW(self):
        return (AxisWindow.cursor_u, AxisWindow.cursor_v, AxisWindow.cursor_w)

    @classmethod
    def FreezeCursorPosition(self, u=-1, v=-1, w=-1):
        """
        Freeze/Unfreeze cursor position (u, v, w)
        """
        if u != -1:
            AxisWindow.is_set_u = u
        if v != -1:
            AxisWindow.is_set_v = v
        if w != -1:
            AxisWindow.is_set_w = w

    @classmethod
    def IsUVWAllSet(self):
        return AxisWindow.is_set_u * AxisWindow.is_set_v * AxisWindow.is_set_w

    @classmethod
    def IsUVWNoneSet(self):
        return ((not AxisWindow.is_set_u) * (not AxisWindow.is_set_v) * (not AxisWindow.is_set_w))


class AxisWindow(MapWindow):
    # Class attributes
    cursor_u = None
    cursor_v = None
    cursor_w = None
    is_set_u = False
    is_set_v = False
    is_set_w = False
    is_size_set = False
    zoom_size = None
    L = None

    def __init__(self, parent, id, axis):
        super(AxisWindow, self).__init__(parent, id)
        self.axis = axis
        # Cursors
        self.we_cursor = wx.StockCursor(wx.CURSOR_SIZEWE)
        self.ns_cursor = wx.StockCursor(wx.CURSOR_SIZENS)

    def RefreshCursor(self):
        if self.axis == "los":
            x = AxisWindow.is_set_u
            y = AxisWindow.is_set_v
            z = AxisWindow.is_set_w
        elif self.axis == "u":
            x = AxisWindow.is_set_w
            y = AxisWindow.is_set_v
            z = AxisWindow.is_set_u
        elif self.axis == "v":
            x = AxisWindow.is_set_u
            y = AxisWindow.is_set_w
            z = AxisWindow.is_set_v

        if (x * y):
            if not z:
                self.SetCursor(self.cross_cursor)
            else:
                if AxisWindow.is_size_set:
                    self.SetCursor(self.arrow_cursor)
                else:
                    self.SetCursor(self.sizing_cursor)
        else:
            if y:
                self.SetCursor(self.we_cursor)
            elif x:
                self.SetCursor(self.ns_cursor)
            else:
                self.SetCursor(self.cross_cursor)

    def MakeBitmap(self):
        """
        Window bitmap creation method.
        """
        # Undefined IMAGE, return
        if self.image is None:
            return

        # Size of the window
        if AxisWindow.L is not None:
            size = AxisWindow.L
        else:
            size = 1.0
        bnx = int(round(self.nx / size))
        bny = int(round(self.ny / size))
        if ((self.img_nx == bnx) * (self.img_ny == bny)):
            image = self.image
        else:
            image = self.image.Scale(bnx, bny)

        self.bitmap = image.ConvertToBitmap()

    def OnPaint(self, event):
        """
        Window painting method.
        """
        # Drawing context
        dc = wx.PaintDC(self)

        # If no defined BITMAP
        if self.bitmap is None:
            # If no IMAGE defined, can't build any bitmap, nothing to do
            if self.image is None:
                return
            # Buil a brand naw bitmap
            else:
                self.MakeBitmap()

        # Cross coordinates (pixels)
        nx = self.nx
        ny = self.ny
        i = None
        iset = False
        j = None
        jset = False
        if self.axis == "los":
            if AxisWindow.cursor_u is not None:
                i = int((AxisWindow.cursor_u + 0.5) * nx)
            if AxisWindow.cursor_v is not None:
                j = int((0.5 - self.cursor_v) * ny)
            iset = AxisWindow.is_set_u
            jset = AxisWindow.is_set_v
        elif self.axis == "u":
            if AxisWindow.cursor_w is not None:
                i = int((0.5 - AxisWindow.cursor_w) * nx)
            if AxisWindow.cursor_v is not None:
                j = int((0.5 - self.cursor_v) * ny)
            iset = AxisWindow.is_set_w
            jset = AxisWindow.is_set_v
        elif self.axis == "v":
            if AxisWindow.cursor_u is not None:
                i = int((AxisWindow.cursor_u + 0.5) * nx)
            if AxisWindow.cursor_w is not None:
                j = int((0.5 - self.cursor_w) * ny)
            iset = AxisWindow.is_set_u
            jset = AxisWindow.is_set_w

        # Background bitmap drawing
        if AxisWindow.L is None:
            dc.DrawBitmap(self.bitmap, 0, 0, useMask=False)
        else:
            istart = i - int(round(i / AxisWindow.L))
            jstart = j - int(round(j / AxisWindow.L))
            dc.DrawBitmap(self.bitmap, istart, jstart, useMask=False)

        if RegionFinderTab.IsUVWAllSet():  # Draw zoom box
            c = AxisWindow.cross_size

            # Draw target little cross
            dc.SetPen(self.redPen)
            dc.DrawLine(i - c, j, i + c, j)
            dc.DrawLine(i, j - c, i, j + c)

            # Draw zoom box
            if AxisWindow.L is None:
                zb = int(0.5 * AxisWindow.zoom_size * nx)
            else:
                zb = int(0.5 * (AxisWindow.zoom_size + 1.0 - AxisWindow.L) * nx)
            if not AxisWindow.is_size_set:
                dc.SetPen(self.whitePen)

            dc.DrawLine(i - zb, j - zb, i - zb, j + zb + 1)
            dc.DrawLine(i - zb, j + zb, i + zb + 1, j + zb)
            dc.DrawLine(i + zb, j - zb, i + zb, j + zb + 1)
            dc.DrawLine(i - zb, j - zb, i + zb + 1, j - zb)
        else:  # Draw full cross
            if i is not None:
                if iset:
                    dc.SetPen(self.redPen)
                else:
                    dc.SetPen(self.whitePen)
                dc.DrawLine(i, 0, i, ny)

            if j is not None:
                if jset:
                    dc.SetPen(self.redPen)
                else:
                    dc.SetPen(self.whitePen)
                dc.DrawLine(0, j, nx, j)
