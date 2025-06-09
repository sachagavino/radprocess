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
#ifndef _ADAPTIVE_GAUSSIAN_BLUR
#define _ADAPTIVE_GAUSSIAN_BLUR
#endif

static void adaptive_gaussian_blur_C(double * map_filtered, double * original_map, int * same_value_pixel_size_map, int map_max_size_i, int map_max_size_j);

static void compute_same_value_pixel_size_map_C(int * result_map, double * original_map, int map_max_size_i, int map_max_size_j);