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
#ifndef _H_RAMSES_UTILS
#define _H_RAMSES_UTILS

#define FELEM2D(A, I, J, NX, NY)          ((A)[(I) + (J)*(NX)])
#define CELEM2D(A, I, J, NX, NY)          ((A)[(I)*(NY) + (J)])
#define CELEM3D(A, I, J, K, NX, NY, NZ)   ((A)[(I)*(NY)*(NZ) + (J)*(NZ) + (K)])
#define FREE_IF_NOTNULL(A)   if(A != NULL) free(A);

#endif
