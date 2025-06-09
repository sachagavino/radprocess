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
#ifndef _H_RAMSES_READ_PARTS
#define _H_RAMSES_READ_PARTS

#include "fio.h"
#include <stdint.h>

#define READ_PARTS_NO_ERROR            0
#define READ_PARTS_MEMORY_ERROR        250


typedef struct RAMSES_PartData_t {
	int ncpu;
	int ndim;
	int npart;
	int nstar;
	int nsink;
	double * pos;  /* (npart, ndim) */
	double * vel;  /* (npart, ndim) */
	double * mass;
	int64_t * id_64;
	int * id_32;
	int * level;
	double * epoch;
	double * metal;
} RAMSES_PartData_t;


RAMSES_PartData_t * RAMSES_PartData_New(void);

void RAMSES_PartData_Free(RAMSES_PartData_t * part);

int RAMSES_PartData_Read(RAMSES_PartData_t * part, double boxlen, FIO_File * file, int long_int);

#endif
