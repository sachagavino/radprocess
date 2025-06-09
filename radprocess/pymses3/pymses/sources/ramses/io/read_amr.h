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
#ifndef _H_RAMSES_READ_AMR
#define _H_RAMSES_READ_AMR

#include "fio.h"

#define MIN(A, B)   (A < B ? A : B)

#define READ_AMR_NO_ERROR             0
#define READ_AMR_MEMORY_ERROR         350
#define RAED_AMR_NGRID_MISMATCH_ERROR 360

typedef struct RAMSES_AmrHeader_t {
	int ndim;                 // Number of dimensions
	int ncpu;                 // Number of CPUs
	int coarse_shape[3];      // Shape of the coarse (level 0) cell block
	int ncoarse;              // Total number of coarse (level 0) cells
	int nboundary;            // Number of boundary regions defined
	int levelmax;             // Maximum AMR level
	int ngridmax;             // Maximum grids for this CPU

	int noutputs;             // Number of predetermined outputs
	double * out_times;       // Sim time for the outputs, if specified
	double * out_aexps;       // Exp. fact. for the outputs, if specified
	double boxlen;            // Size of the simulation box

	int iout;                 // Current output ID
	double aexp;              // Current expansion factor
	double time;              // Current sim time
	int ngrids_lmax;
} RAMSES_AmrHeader_t;

typedef struct RAMSES_AmrStruct_t {
	int ndim;
	int twotondim;            // (1 << ndim)
	int ncoarse;
	int ncpu;
	int nboundary;
	int levelmax;
	int ngrids;
	int ngrids_lmax;
	int readlmax;

	/* Allocated grids, by region and level.
	 *     shape: (ncpu+nboundary, levelmax)
	 */
	int * ngridlevel;
	int * ngridbound;

	int * grid_indices;

	/* Coarse level arrays */
	/* Sons of the coarse cells
	 *     shape: (ncoarse,)
	 */
	int * coarse_son_indices;

	/* Grid-based arrays */

	/* Positions of grid centers.
	 *     shape: (ngrids, ndim);
	 */
	double * grid_centers;

	/* Cell-based arrays */

	/* Indices of the son grids
	 *     shape: (ngrids, twotondim)
	 */
	int * son_indices;

} RAMSES_AmrStruct_t;

/* AmrHeader */
RAMSES_AmrHeader_t * RAMSES_AmrHeader_New(void);
void RAMSES_AmrHeader_Free(RAMSES_AmrHeader_t * header);
int RAMSES_AmrHeader_Read(RAMSES_AmrHeader_t * header, FIO_File * file);


/* AmrStruct */
RAMSES_AmrStruct_t * RAMSES_AmrStruct_New(void);
void RAMSES_AmrStruct_Free(RAMSES_AmrStruct_t * amr);


/* General AMR functions */
int RAMSES_ReadAmr(RAMSES_AmrHeader_t * header, RAMSES_AmrStruct_t * amr, const char * amr_filename, int readlmax,
                   int swap);

#endif
