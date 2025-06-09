/******************************************************************************************
*   This file is part of PyMSES.                                                          *
*                                                                                         *
*   PyMSES is free software: you can redistribute it and/or modify                        *
*   it under the terms of the GNU General Public License as published by                  *
*   the Free Software Foundation, either version 3 of the License, or                     *
*   (at your option) any later version.                                                   *
*                                                                                         *
*   PyMSES is distributed in the hope that it will be useful,                             *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of                        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                         *
*   GNU General Public License for more details.                                          *
*                                                                                         *
*   You should have received a copy of the GNU General Public License                     *
*   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.                       *
******************************************************************************************/
#ifndef _RAY_TRACE_C_CODE_H
#define _RAY_TRACE_C_CODE_H
#endif

#include<math.h> 

static void ray_trace_full_octree(double * I, double * ray_origins, double * ray_vects,
	double * grid_centers, double * cell_centers, int * sons, int * cell_levels,
	double * scalar_data, int * neighbors, double * ray_lengths,
	long * param_info_int, double gamma, double * vlines, double * wlines,
	double * rlines, double * glines, double * blines, double * alines);
