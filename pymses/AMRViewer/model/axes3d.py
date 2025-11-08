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

import numpy as N


class Axes3DModel(object):
    axis = {"x": N.array([1.0, 0.0, 0.0]),
            "-x": N.array([-1.0, 0.0, 0.0]),
            "y": N.array([0.0, 1.0, 0.0]),
            "-y": N.array([0.0, -1.0, 0.0]),
            "z": N.array([0.0, 0.0, 1.0]),
            "-z": N.array([0.0, 0.0, -1.0])}

    def __init__(self):
        """
        Axes3d : creates a set of points (3D coordinates) to draw the x, y and z axis
        on the view.
        """
        # points
        self.points = []
        df = 0.1
        # x-arrow
        self.points.append(Axes3DModel.axis["x"])
        self.points.append(N.array([1.0 - df, df / 2, 0.0]))
        self.points.append(N.array([1.0 - df, -df / 2, 0.0]))
        # y-arrow
        self.points.append(Axes3DModel.axis["y"])
        self.points.append(N.array([0.0, 1.0 - df, -df / 2]))
        self.points.append(N.array([0.0, 1.0 - df, df / 2]))
        # z-arrow
        self.points.append(Axes3DModel.axis["z"])
        self.points.append(N.array([df / 2, 0.0, 1.0 - df]))
        self.points.append(N.array([-df / 2, 0.0, 1.0 - df]))

    def RefreshPoints(self):
        """
        Refresh the (u,v) coordinates of the point list according
        to the current self.uvect and self.vvect and stores it into
        self.uvpoints
        """
        n = len(self.points)
        self.uvpoints = []
        for i in range(n):
            u = N.sum(self.points[i] * self.uvect)
            v = N.sum(self.points[i] * self.vvect)
            self.uvpoints.append(N.array([u, v]))

    def GetUVPoints(self):
        return self.uvpoints

    def UpdateUVVectors(self, rot_2D):
        if ((self.du == 0.0) * (self.dv == 0.0)):
            return
        # Rotation distance
        I_move = N.sqrt(self.du ** 2 + self.dv ** 2)
        # Rotation axis
        wvect = self.GetLosAxis()
        if rot_2D:  # We rotate in the current u-v plane
            uc = self.ucursor
            vc = self.vcursor
            uo = uc - self.du
            vo = vc - self.dv
            normc = N.sqrt(uc ** 2 + vc ** 2)
            normo = N.sqrt(uo ** 2 + vo ** 2)
            uo = uo / normo
            vo = vo / normo
            uc = uc / normc
            vc = vc / normc
            stheta = uo * vc - uc * vo
            ctheta = uo * uc + vo * vc
            # Rotation angle
            alpha_rad = -N.arctan2(stheta, ctheta)
            move_axis = wvect
        else:  # free 3D rotation
            # Rotation angle
            alpha_rad = I_move

            move_vect = -(self.du * self.uvect + self.dv * self.vvect) / I_move
            move_vect = move_vect / N.linalg.norm(move_vect, 2)
            move_axis = N.cross(wvect, move_vect)

        marr = move_axis
        # P = rot_axis * transpose(rot_axis) => Projection along the
        # rotation axis
        self.Prot[:, 0] = marr * marr[0]
        self.Prot[:, 1] = marr * marr[1]
        self.Prot[:, 2] = marr * marr[2]

        # Q => antisymetrical representation of the rotation axis
        self.Qrot[:, :] = 0.0
        self.Qrot[0, 1] = -marr[2]
        self.Qrot[1, 0] = marr[2]
        self.Qrot[0, 2] = marr[1]
        self.Qrot[2, 0] = -marr[1]
        self.Qrot[1, 2] = -marr[0]
        self.Qrot[2, 1] = marr[0]

        ma = self.Prot + (self.Irot - self.Prot) * N.cos(alpha_rad) + self.Qrot * N.sin(alpha_rad)

        # Transform the u and v vectors
        au = self.uvect
        av = self.vvect
        au2 = N.tensordot(au, ma, axes=[0, 1])
        self.uvect = au2 / N.linalg.norm(au2, ord=2)
        av2 = N.tensordot(av, ma, axes=[0, 1])
        av2 = av2 - N.dot(av2, self.uvect) * self.uvect
        self.vvect = av2 / N.linalg.norm(av2, ord=2)

        self.du = 0.
        self.dv = 0.

    def MoveCursor(self, dx, dy, newx, newy):
        if ((dx != 0.0) + (dy != 0)):
            self.ucursor = newx
            self.vcursor = newy
            self.du = self.du + dx
            self.dv = self.dv + dy

    def GetLosAxis(self):
        """
        Returns the coordinates of the line-of-sight axis.
        """
        wvect = N.cross(self.uvect, self.vvect)
        wvect = wvect / N.linalg.norm(wvect, 2)
        return wvect

    def GetUVAxes(self):
        return (self.uvect, self.vvect)

    def SetUVAxis(self, u, v):
        if u in list(Axes3DModel.axis.keys()):
            self.uvect = Axes3DModel.axis[u]
        else:
            self.uvect = u
        if v in list(Axes3DModel.axis.keys()):
            self.vvect = Axes3DModel.axis[v]
        else:
            self.vvect = v

    def reset(self):
        print("Reset Axes3DModel")
        self.du = 0.0
        self.dv = 0.0

        # Sets default u and v axis
        self.SetUVAxis(u="x", v="y")
        self.RefreshPoints()

        self.Qrot = N.zeros(shape=(3, 3))
        self.Prot = N.zeros(shape=(3, 3))
        self.Irot = N.eye(3)
