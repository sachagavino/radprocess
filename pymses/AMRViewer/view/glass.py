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


class MagGlass(wx.Window):
    mag_factor = 8

    def __init__(self, parent):
        WID = wid()
        self.redPen = wx.Pen("RED")
        super(MagGlass, self).__init__(parent, WID.GLASS_WIN, style=wx.RAISED_BORDER, size=(200, 200))
        self.SetBackgroundColour(wx.BLACK)

        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.bitmap = None

    def OnPaint(self, event):
        if self.bitmap is None:
            return
        dc = wx.PaintDC(self)
        dc.DrawBitmap(self.bitmap, 0, 0, useMask=False)

        i = 100
        j = 100
        dc.SetPen(self.redPen)
        mgf = MagGlass.mag_factor
        # draw a thick re-sized cross to show selected pixel value
        shift = 0
        if mgf > 4:
            shift = 1
        if mgf % 2 != 0 or mgf == 2:
            # odd number of pixel : add more pixels
            for d in range(mgf / 4, mgf * 3 / 4 + 1):
                dc.DrawLine(i - mgf, j + d + shift, i + 2 * mgf, j + d + shift)
                dc.DrawLine(i + d, j - mgf + shift, i + d, j + 2 * mgf + shift)
        else:
            for d in range(mgf / 4, mgf * 3 / 4):
                dc.DrawLine(i - mgf, j + d + shift, i + 2 * mgf, j + d + shift)
                dc.DrawLine(i + d, j - mgf + shift, i + d, j + 2 * mgf + shift)
