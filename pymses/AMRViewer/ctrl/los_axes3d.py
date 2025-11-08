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
import wx


def get_unit_coords_from_pixel_coords(x, y, nx, ny):
    a = 2. * float(x) / nx - 1.0
    b = 1.0 - 2. * float(y) / ny
    return (a, b)


class LineOfSightController(object):
    timing_ms = 20  # Refreshing time interval

    def __init__(self):
        self.a3d_time = 0

    def OnShortcutAxis(self, event):
        WID = wid()
        shortcut_id = event.GetId()
        if shortcut_id == WID.LOS_TOOL_X:
            self.model.SetUVAxis(u="y", v="z")
        elif shortcut_id == WID.LOS_TOOL_MX:
            self.model.SetUVAxis(u="z", v="y")
        elif shortcut_id == WID.LOS_TOOL_Y:
            self.model.SetUVAxis(u="z", v="x")
        elif shortcut_id == WID.LOS_TOOL_MY:
            self.model.SetUVAxis(u="x", v="z")
        elif shortcut_id == WID.LOS_TOOL_Z:
            self.model.SetUVAxis(u="x", v="y")
        elif shortcut_id == WID.LOS_TOOL_MZ:
            self.model.SetUVAxis(u="y", v="x")
        else:
            uvect, vvect = self.model.GetUVAxes()
            wvect = self.model.GetLosAxis()
            if shortcut_id == WID.LOS_TOOL_ROT_LEFT:
                self.model.SetUVAxis(u=-wvect, v=vvect)
            elif shortcut_id == WID.LOS_TOOL_ROT_RIGHT:
                self.model.SetUVAxis(u=wvect, v=vvect)
            elif shortcut_id == WID.LOS_TOOL_ROT_TOP:
                self.model.SetUVAxis(u=uvect, v=wvect)
            elif shortcut_id == WID.LOS_TOOL_ROT_BOTTOM:
                self.model.SetUVAxis(u=uvect, v=-wvect)

        self.axes3D_refresh()
        if self.A3d_autorefresh:
            self.UpdateImages()

    def A3dOnMouse(self, event):
        a3d = event.GetEventObject()
        refresh = False
        if event.Entering():
            refresh = a3d.SetMouseIn(True)
        elif event.Leaving():
            refresh = a3d.SetMouseIn(False)
        elif event.LeftDown():
            a3d.SetClickDown(True)
            a3d.UpdateCursor()
            nx, ny = a3d.GetSize()
            x, y = event.GetPosition()
            cursor_x, cursor_y = get_unit_coords_from_pixel_coords(x, y, nx, ny)
            a3d.SetCursorPosition(cursor_x, cursor_y)
        elif event.LeftUp():
            a3d.SetClickDown(False)
            a3d.Refresh()
            a3d.UpdateCursor()
            if self.A3d_autorefresh:
                self.UpdateImages()
        elif event.Dragging():
            if event.LeftIsDown():
                nx, ny = a3d.GetSize()
                i, j = event.GetPosition()
                x, y = get_unit_coords_from_pixel_coords(i, j, nx, ny)
                x0, y0 = a3d.GetPos()
                a3d.SetCursorPosition(x, y)
                d_x = x - x0
                d_y = y - y0
                self.model.MoveCursor(d_x, d_y, x, y)
                t = event.GetTimestamp()
                if (t - self.a3d_time) > LineOfSightController.timing_ms:
                    self.a3d_time = t
                    self.model.UpdateUVVectors(rot_2D=event.ControlDown())
                    refresh = True
        if refresh:
            self.axes3D_refresh()

    def A3dOnKeyPress(self, event):
        WID = wid()
        k = event.GetKeyCode()
        if k == wx.WXK_CONTROL:
            a3d = self.widgets[WID.EXPANDER_LOS].GetAxes3d()
            a3d.SetCtrlDown(True)
            a3d.UpdateCursor()
        event.Skip()

    def A3dOnKeyRelease(self, event):
        WID = wid()
        k = event.GetKeyCode()
        if k == wx.WXK_CONTROL:
            a3d = self.widgets[WID.EXPANDER_LOS].GetAxes3d()
            a3d.SetCtrlDown(False)
            a3d.UpdateCursor()
        event.Skip()

    def axes3D_refresh(self):
        WID = wid()
        self.model.RefreshPoints()
        self.widgets[WID.EXPANDER_LOS].RefreshAxes3D()
