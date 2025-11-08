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


class RegionFinderController(object):
    timing_ms = 20  # Refreshing time interval

    def __init__(self):
        self.rf_time = 0

    def RegionFinderOnMouse(self, event):
        """
        Capture and handles mouse events occuring in the RegionFinder tab windows (line-of-sight, u and -v views).
        """
        WID = wid()
        # RegionFinderTab view
        rf_view = self.widgets[WID.RF_TAB]

        # Axis window and its ID/axis name
        win = event.GetEventObject()
        win_id = win.GetId()
        axis = win.axis

        # nothing to do if the window doesn't have a map to display
        if win.bitmap is None:
            return

        t = event.GetTimestamp()
        if ((event.Moving() + event.Dragging()) * (
                    (t - self.rf_time) < RegionFinderController.timing_ms)):  # Refreshment frequency
            return
        self.rf_time = t

        make_bitmap = False
        if event.Leaving():
            # Leaving the window
            if axis == "los":
                rf_view.MoveCursor(u=None, v=None)
            elif axis == "u":
                rf_view.MoveCursor(v=None, w=None)
            elif axis == "v":
                rf_view.MoveCursor(u=None, w=None)
        else:
            if event.RightUp():
                # Right click : UNSET behavior
                if rf_view.IsZoomSizeIsSet():  # Unset zoom size
                    rf_view.SetSizeSet(False)

                    # Refresh cursor types on all windows
                    rf_view.RefreshWindowCursors()
                else:  # Unfreeze cursor position
                    if not rf_view.IsUVWNoneSet():
                        rf_view.ResetZoomBox()
                        make_bitmap = True
                        if axis == "los":
                            rf_view.FreezeCursorPosition(u=False, v=False)
                        elif axis == "u":
                            rf_view.FreezeCursorPosition(v=False, w=False)
                        elif axis == "v":
                            rf_view.FreezeCursorPosition(u=False, w=False)
                        # Refresh cursor types on all windows
                        rf_view.RefreshWindowCursors()

                        self.model.SetZoomSize(rf_view.GetZoomSize())

            # Cursor movement
            # Pointer position of the event :
            #  - (i, j) = pixel bitmap coordinates
            #  - (x, y) = [-0.5, 0.5]^2 coordinates ([0.5, 0] pointing right and [0,0.5] pointing upward)
            i, j = event.GetPosition()
            x = float(i) / win.nx - 0.5
            y = 0.5 - float(j) / win.ny
            # Move cursor position
            if axis == "los":
                rf_view.MoveCursor(u=x, v=y)
            elif axis == "u":
                rf_view.MoveCursor(v=y, w=-x)
            elif axis == "v":
                rf_view.MoveCursor(u=x, w=y)

            if event.LeftUp():
                # Left click : SET behavior
                if not rf_view.IsUVWAllSet():  # Freeze (previously unfreezed) cursor position
                    if axis == "los":
                        rf_view.FreezeCursorPosition(u=True, v=True)
                    elif axis == "u":
                        rf_view.FreezeCursorPosition(v=True, w=True)
                    elif axis == "v":
                        rf_view.FreezeCursorPosition(u=True, w=True)

                    if rf_view.IsUVWAllSet():  # All positions are freezed, init. zoom box
                        rf_view.InitZoomSize()
                        self.model.SetZoomSize(rf_view.GetZoomSize())
                else:  # All positions are already freezed => freeze zoom box
                    if not rf_view.IsZoomSizeIsSet():
                        rf_view.SetSizeSet(True)
                    else:
                        self.UpdateImages()

                # Refresh cursor types on all windows
                rf_view.RefreshWindowCursors()

            # Mouse wheel behavior to modify zoom box size
            wh = event.GetWheelRotation()
            # Wheel event to modify the region size
            if (wh != 0) * rf_view.IsUVWAllSet() * (not rf_view.IsZoomSizeIsSet()):
                make_bitmap = rf_view.ChangeZoomSize(wh)
                self.model.SetZoomSize(rf_view.GetZoomSize())

        u, v, w = rf_view.GetCursorUVW()
        self.model.UpdateCursorUVWPosition(u, v, w)
        rf_view.RefreshWholeTab(make_bitmap)

    def RegionFinderButtonClick(self, do_all_processing, verbose):
        WID = wid()
        print("Start updating view ...")
        # RegionFinderTab view
        rf_view = self.widgets[WID.RF_TAB]
        if self.model.ro is not None:
            nxy_list = rf_view.GetWindowSizesDict()
            proceed = True
            if do_all_processing:
                proceed = self.model.UpdateRegionFinderMaps(nxy_list)
            if proceed:
                self.model.UpdateImage(verbose=verbose)
                rf_view.reset_tab()
                rf_view.UpdateRegionFinderImages()

                # Refresh cursor types on all windows
                rf_view.RefreshWindowCursors()
                self.model.SetZoomSize(rf_view.GetZoomSize())
