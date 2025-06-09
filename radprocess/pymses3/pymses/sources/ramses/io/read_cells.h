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
#ifndef _H_RAMSES_READ_CELLS
#define _H_RAMSES_READ_CELLS

#include "fio.h"
#include "read_amr.h"

#define READ_CELLS_NO_ERROR            0
#define READ_CELLS_MEMORY_ERROR        200
#define READ_CELLS_NODATA_ERROR        201
#define READ_CELLS_NDIM_INVALID        202
#define READ_CELLS_INVALID_IDATA_ERROR 203
#define READ_CELLS_READ_NGRID_MISS     204
#define READ_CELLS_READ_NGRID_OVERFLOW 205
#define READ_CELLS_UNAVAIL_IDATA_ERROR 300

/* Cell data structure definition */
typedef struct RAMSES_CellData_t {
	double * data;
	int nvar_file;
	int nvar_data;
	int * ivar_arr;
	int ncpu;
	int ndim;
	int ngrids;
	int nboundary;
	int levelmax;
} RAMSES_CellData_t;

/* Cell data hydro file header structure definition */
typedef struct RAMSES_CellHydroHeader_t {
	int nvar_file;
	int ncpu;
	int ndim;
	int levelmax;
	double gamma;
} RAMSES_CellHydroHeader_t;

/* Method definitions */
RAMSES_CellData_t * RAMSES_CellData_New(int nvarout, int ndim, int ngrids);
void RAMSES_CellData_Free(RAMSES_CellData_t * cell);
int RAMSES_CellData_Read(RAMSES_CellData_t * cell,	const char * amrdata_filename, int is_hydro, int readlmax, int swap,
                         int grav_compat);

RAMSES_CellHydroHeader_t * RAMSES_CellHydro_Header_New(void);
void RAMSES_CellHydro_Header_Free(RAMSES_CellHydroHeader_t * hdr);
int RAMSES_CellHydro_Header_Read(RAMSES_CellHydroHeader_t * header,	const char * hydrodata_filename, int swap);

#endif

