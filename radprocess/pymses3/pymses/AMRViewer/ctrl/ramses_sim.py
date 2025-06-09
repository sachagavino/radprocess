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


class RamsesSimController(object):
    def __init__(self):
        self.ramses_look_dir = "./"
        self.WID = wid()

    def SelectRamsesDir(self, out_dir):
        """
        Send a SearchValidOutputs request to the model
        and return the returned value
        """
        ok, selected_iout_index = self.model.SearchValidOutputs(out_dir)
        if ok:
            self.widgets[self.WID.RF_TAB].RefreshOutputInfo()
            self.SelectRamsesOutput(selected_iout_index)
        return ok

    def OnSelectIout(self, event):
        iout = int(self.widgets[self.WID.RF_TAB].iout_cb.GetValue())
        try:
            iout_index = self.model.iout_list.index(iout)
        except Exception:
            print("Error : Output ", iout, " is not available !")
            iout_index = self.model.selected_iout_index
            self.widgets[self.WID.RF_TAB].SetRamsesOutput()
        if iout_index != self.model.selected_iout_index:
            self.SelectRamsesOutput(iout_index)

    def RefreshRamsesIoutList(self, event):
        """
        Refresh the ioutput_list.
        - Send a SearchValidOutputs request to the RAMSES model (with a
        None output directory argument).
        - If the result is True,
        """
        if self.model.SearchValidOutputs(None)[0]:
            self.widgets[self.WID.RF_TAB].RefreshOutputInfo(select_last_iout=True)

    def SelectRamsesOutput(self, iout_index=None, iout=None):
        """
        Selects a RAMSES output within the list of available output of the
        current simulation directory.

        Sends a request to the model to select the given output and makes the view refresh itself afterwards.
        """
        if self.model.SetRamsesOutput(iout_index=iout_index, iout=iout):
            self.widgets[self.WID.RF_TAB].SetRamsesOutput()
            self.model.map_engine_dict.clear()
        self.UpdateLogScale()
        self.UpdateFraction()
