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
from .utils import get_icon
import wx.lib.buttons as buttons
from numpy import sqrt
from .dialogs import CameraPNGFileDialog


class MapWindow(wx.Window):
    cross_size = 4

    def __init__(self, parent, id=None):
        WID = wid()
        if id is None:
            id = WID.ANY
        super(MapWindow, self).__init__(parent, id, style=wx.NO_BORDER, size=(255, 255))
        self.SetBackgroundColour(wx.BLACK)
        self.nx, self.ny = self.GetSize()

        # Pens
        self.whitePen = wx.Pen("WHITE")
        self.redPen = wx.Pen("RED")
        # Cursors
        self.no_cursor = wx.StockCursor(wx.CURSOR_BLANK)
        self.cross_cursor = wx.StockCursor(wx.CURSOR_CROSS)
        self.arrow_cursor = wx.StockCursor(wx.CURSOR_ARROW)
        self.sizing_cursor = wx.StockCursor(wx.CURSOR_SIZING)

        self.reset()

        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_SIZE, self.OnResize)

    def reset(self):
        self.SetBackgroundImage(None)
        self.bitmap = None
        self.SetCursor(self.arrow_cursor)

    def SetBackgroundImage(self, pil_image):
        """
        Change Window background image.

        pil_image  ---  new PIL image to set. If None, reset image, nx, ny and bitmpap to None
        """
        if pil_image is None:
            self.image = None
            self.pil_image = None
            self.img_nx = None
            self.img_ny = None
            self.bitmap = None
        else:
            self.pil_image = pil_image
            self.img_nx, self.img_ny = pil_image.size
            self.image = wx.EmptyImage(self.img_nx, self.img_ny)
            self.image.SetData(pil_image.convert("RGB").tostring())
            self.image.SetAlphaData(self.pil_image.tostring()[3::4])

    def RefreshWindow(self, make_bitmap=False):
        if make_bitmap:
            self.MakeBitmap()
        self.Refresh()

    def RefreshCursor(self):
        raise NotImplementedError()

    def MakeBitmap(self):
        """
        Window bitmap creation method.
        """
        # Undefined IMAGE, return
        if self.image is None:
            return

        simg = self.GetScaledImage()
        self.bitmap = simg.ConvertToBitmap()

    def GetScaledImage(self, pil=False):
        """
        Returns the scaled image (to fit the size of the window)
        """
        if self.image is None:
            return None
        # Size of the window
        bnx = self.nx
        bny = self.ny
        if ((self.img_nx == bnx) * (self.img_ny == bny)):
            if pil:
                return self.pil_image
            else:
                return self.image
        else:
            if pil:
                return self.pil_image.resize((bnx, bny))
            else:
                return self.image.Scale(bnx, bny)

    def OnPaint(self, event):
        raise NotImplementedError()

    def OnResize(self, event):
        self.nx, self.ny = self.GetSize()
        self.RefreshWindow(make_bitmap=True)


class MapTab(wx.Panel):
    def __init__(self, parent, name, model, controller):
        WID = wid()
        id = WID.TAB_IDS[name]
        super(MapTab, self).__init__(parent, id)
        self.name = name
        self.model = model

        page_sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(page_sizer)
        b = wx.StaticBox(self, WID.ANY, "line-of-sight view")
        bs = wx.StaticBoxSizer(b, wx.HORIZONTAL)

        # Map window

        self.window = GenericTabMapWindow(self, name)
        self.window.Bind(wx.EVT_MOUSE_EVENTS, controller.TabMapWindowOnMouse)
        bs.Add(self.window, 1, wx.EXPAND | wx.SHAPED | wx.ALIGN_CENTER | wx.ALL, 1)

        page_sizer.Add(bs, 1, wx.EXPAND | wx.ALL, 4)

    def reset(self):
        self.window.reset()
        self.RefreshWholeTab()

    def GetWindowSize(self):
        return self.window.GetSize()

    def UpdateTabImage(self, verbose=True):
        img = self.model.GetTabImage(self.name, verbose=verbose)
        self.window.SetBackgroundImage(img)

    def RefreshWholeTab(self, make_bitmap=False):
        """
        Refresh the tab content
        """
        self.window.RefreshWindow(make_bitmap)


class GenericTabMapWindow(MapWindow):
    mag_box_size = 25  # Size of the magnifier

    def __init__(self, parent, name):
        WID = wid()
        super(GenericTabMapWindow, self).__init__(parent, id=WID.TAB_WIN_IDS[name])
        self.name = name

    def reset(self):
        MapWindow.reset(self)
        self.old_cursor_i = None
        self.old_cursor_j = None
        self.cursor_i = None
        self.cursor_j = None

    def GetMagnifierBoxImage(self):
        """
        Returns the (cropped) image corresponding to the content of the magnifier box
        """
        if ((self.cursor_i is None) + (self.cursor_j is None)):
            return None
        spimg = self.GetScaledImage(pil=True)
        nb = GenericTabMapWindow.mag_box_size
        img = spimg.crop(
            (self.cursor_i - nb / 2, self.cursor_j - nb / 2, self.cursor_i + nb / 2, self.cursor_j + nb / 2))
        return img

    def GetRuleDist(self):
        if ((self.old_cursor_i is not None) * (self.old_cursor_j is not None) * \
                    (self.cursor_i is not None) * (self.cursor_j is not None)):
            di2 = (self.cursor_i - self.old_cursor_i) ** 2
            dj2 = (self.cursor_j - self.old_cursor_j) ** 2
            dist = sqrt(float(di2 + dj2)) / self.nx
            return dist
        else:
            return None

    def RefreshCursor(self):
        if self.bitmap is not None:
            self.SetCursor(self.cross_cursor)
        else:
            self.SetCursor(self.arrow_cursor)

    def MoveCursor(self, i, j):
        self.cursor_i = i
        self.cursor_j = j

    def FreezeOrigin(self, i, j):
        self.old_cursor_i = i
        self.old_cursor_j = j

    def GetMapCursorPosition(self):
        """
        Return the indices of the map cell under the cursor,
        (None, None) if the cursor is out of the window
        """
        if ((self.cursor_i is None) + (self.cursor_j is None)):
            return (None, None)
        else:
            if ((self.nx == self.img_nx) * (self.ny == self.img_ny)):
                return (self.cursor_i, self.cursor_j)
            else:
                x = (self.cursor_i / float(self.nx) * float(self.img_nx))
                y = (self.cursor_j / float(self.ny) * float(self.img_ny))
                return (x, y)

    def OnPaint(self, event):
        """
        Window painting method.
        """
        # Drawing context
        dc = wx.PaintDC(self)

        # If no defined BITMAP
        if self.bitmap is None:
            # If no IMAGE defined, can't buil any bitmap, nothing to do
            if self.image is None:
                return
            # Buil a brand new bitmap
            else:
                self.MakeBitmap()

        # Background bitmap drawing
        dc.DrawBitmap(self.bitmap, 0, 0, useMask=False)

        # Cross coordinates (pixels)
        nx = self.nx
        ny = self.ny
        i = self.cursor_i
        j = self.cursor_j
        oi = self.old_cursor_i
        oj = self.old_cursor_j

        # Draw target little cross
        if ((i is not None) * (j is not None)):
            # rule line:
            dc.SetPen(self.redPen)
            c = MapWindow.cross_size
            dc.DrawLine(i - c, j, i + c + 1, j)
            dc.DrawLine(i, j - c, i, j + c + 1)
            if ((oi is not None) * (oj is not None)):
                dc.DrawLine(oi - c, oj, oi + c + 1, oj)
                dc.DrawLine(oi, oj - c, oi, oj + c + 1)
                dc.SetPen(self.whitePen)
                dc.DrawLine(i, j, oi, oj)

            # magnifier frame:
            l = GenericTabMapWindow.mag_box_size / 2
            dc.SetPen(self.whitePen)
            dc.DrawLine(i - l, j - l, i + l + 1, j - l)
            dc.DrawLine(i + l, j - l, i + l, j + l + 1)
            dc.DrawLine(i - l, j + l, i + l + 1, j + l)
            dc.DrawLine(i - l, j - l, i - l, j + l + 1)
            # dc.SetPen(self.redPen)
            # dc.DrawLine(i-1,j-1, i+2, j-1)
            # dc.DrawLine(i+1,j-1, i+1, j+2)
            # dc.DrawLine(i-1,j+1, i+2, j+1)
            # dc.DrawLine(i-1,j-1, i-1, j+2)

            # TODO
