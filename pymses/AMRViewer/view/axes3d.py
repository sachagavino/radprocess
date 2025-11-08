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
from math import sqrt
from .utils import get_icon


class Axes3D(wx.Window):
    ctrl_down = False
    click_down = False
    cursor_x = None
    cursor_y = None
    is_mouse_in = False

    def __init__(self, parent):
        """
        Axes3D window where the u, v and x, y, z axis are displayed.
        The user has a drag acces to these axis to define his own line-of-sight axis interactively.

        parent --- parent object of the Axes3D window.
        """
        WID = wid()
        super(Axes3D, self).__init__(parent, WID.LOS_A3D, style=wx.WANTS_CHARS, size=(200, 200))
        self.SetBackgroundColour(wx.WHITE)
        self.Bind(wx.EVT_PAINT, self.OnPaint)

        # Cursors
        img = wx.Image(get_icon("rotate_2D.png"), type=wx.BITMAP_TYPE_PNG)
        img.ConvertAlphaToMask()
        img.SetOptionInt(wx.IMAGE_OPTION_CUR_HOTSPOT_X, 10)
        img.SetOptionInt(wx.IMAGE_OPTION_CUR_HOTSPOT_Y, 10)
        self.rot_2D_cursor = wx.CursorFromImage(img)
        self.cursor_sizing = wx.StockCursor(wx.CURSOR_SIZING)
        self.cursor_hand = wx.StockCursor(wx.CURSOR_HAND)

        # Pens
        self.redPen = wx.Pen("RED")
        self.blackPen = wx.Pen("BLACK")
        self.bluePen = wx.Pen('#0099f7')
        self.preview_image = None

    def SetMouseIn(self, is_mouse_in):
        Axes3D.is_mouse_in = is_mouse_in
        if is_mouse_in:
            self.SetFocus()
        return True

    def SetCursorPosition(self, x, y):
        Axes3D.cursor_x = x
        Axes3D.cursor_y = y

    def SetClickDown(self, is_clicked):
        Axes3D.click_down = is_clicked

    def SetCtrlDown(self, is_ctrl_down):
        Axes3D.ctrl_down = is_ctrl_down

    def GetPos(self):
        return (Axes3D.cursor_x, Axes3D.cursor_y)

    def UpdateCursor(self):
        if Axes3D.click_down:
            if Axes3D.ctrl_down:
                self.SetCursor(self.rot_2D_cursor)
            else:
                self.SetCursor(self.cursor_sizing)
        else:
            self.SetCursor(self.cursor_hand)

    def SetUVPoints(self, points):
        self.point_list = points

    def SetPreviewImage(self, preview_image_bitmap):
        self.preview_image = preview_image_bitmap

    def OnPaint(self, event):
        """
        Axes3D view painting method
        """
        # Window size
        nx, ny = self.GetSize()
        # Window drawing context
        dc = wx.PaintDC(self)

        # Paint border when mouse is over the Window
        if (Axes3D.is_mouse_in + Axes3D.click_down):
            dc.SetPen(self.bluePen)
            dc.DrawRectangle(0, 0, nx, ny)

        # draw preview cube
        if self.preview_image is not None:
            dc.DrawBitmap(self.preview_image, 0, 0)

        # Paint axes
        # Center of the window coordinates
        ix = int(nx / 2.)
        iy = int(ny / 2.)
        # Pen
        dc.SetTextForeground(wx.BLACK)
        dc.SetPen(self.blackPen)
        # Shrinking factor
        sf = 0.8
        # Axis name list
        tlist = ["x", "y", "z"]
        # Gap between tip of arrow and axis name text pixels)
        tgap = 3
        for i in range(3):
            # Draw axis line
            coords = self.point_list[3 * i]
            ix2 = int((sf * coords[0] + 1.) / 2. * nx)
            iy2 = int((1. - sf * coords[1]) / 2. * ny)
            r2 = sqrt(float((ix2 - ix) ** 2 + (iy2 - iy) ** 2))
            dc.DrawLine(ix, iy, ix2, iy2)

            # Draw arrow lines
            for j in range(1, 3):
                coords2 = self.point_list[3 * i + j]
                ix3 = int((sf * coords2[0] + 1.) / 2. * nx)
                iy3 = int((1. - sf * coords2[1]) / 2. * ny)
                dc.DrawLine(ix2, iy2, ix3, iy3)

            # Axis name positionning/drawing
            t = tlist[i]
            # Axis name text extent (pixels)
            t_nx, t_ny = dc.GetTextExtent(t)
            # Axis name text radius (pixels)
            rt = sqrt((t_nx ** 2 + t_ny ** 2) / 4.)
            # Axis name text position
            if r2 < 1.:
                ixt = int(ix - (tgap + rt) - (t_nx / 2.))
                iyt = int(iy + (tgap + rt) - (t_ny / 2.))
            else:
                ixt = int(ix + (ix2 - ix) * ((r2 + tgap + rt) / r2) - (t_nx / 2.))
                iyt = int(iy + (iy2 - iy) * ((r2 + tgap + rt) / r2) - (t_ny / 2.))
            # Draw axis name text
            dc.DrawTextPoint(t, (ixt, iyt))

        # U and V axis
        dc.SetTextForeground(wx.RED)
        dc.SetPen(self.redPen)
        ix0 = int(0.05 * nx)
        iy0 = int(0.95 * ny)
        # U axis
        t = "u"
        t_nx, t_ny = dc.GetTextExtent(t)
        ixu = ix0 + int(0.15 * nx)
        iyu = iy0
        dc.DrawTextPoint(t, (ixu + tgap, iyu - int(t_ny / 2.)))
        dc.DrawLine(ix0, iy0, ixu, iyu)
        dc.DrawLine(ixu, iyu, ixu - 5, iyu - 2)
        dc.DrawLine(ixu, iyu, ixu - 5, iyu + 2)
        # V axis
        t = "v"
        t_nx, t_ny = dc.GetTextExtent(t)
        ixv = ix0
        iyv = iy0 - int(0.15 * nx)
        dc.DrawTextPoint(t, (ixv - int(t_nx / 2.), iyv - tgap - t_ny))
        dc.DrawLine(ix0, iy0, ixv, iyv)
        dc.DrawLine(ixv, iyv, ixv - 2, iyv + 5)
        dc.DrawLine(ixv, iyv, ixv + 2, iyv + 5)
