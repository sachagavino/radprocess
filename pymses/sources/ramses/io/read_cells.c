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
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "fio.h"
#include "utils.h"
#include "read_cells.h"

RAMSES_CellData_t * RAMSES_CellData_New(int nvarout, int ndim, int ngrids)
{
	RAMSES_CellData_t * cell;
	cell = malloc(sizeof(RAMSES_CellData_t));
	if (cell != NULL)
	{
	    cell->data = NULL;
	    cell->nvar_file = 0;
		cell->nvar_data = nvarout;
		cell->ivar_arr = malloc( ((int)nvarout) * sizeof(int) );
		if (cell->ivar_arr == NULL)
		{
			RAMSES_CellData_Free(cell);
			return NULL;
		}
	    cell->ncpu = 0;
	    cell->ndim = ndim;
		cell->ngrids = ngrids;
		cell->nboundary = 0;
		cell->levelmax = 0;
	}
	return cell;
}


void RAMSES_CellData_Free(RAMSES_CellData_t * cell)
{
	if (cell != NULL)
	{
		FREE_IF_NOTNULL(cell->ivar_arr);
		free(cell);
	}
}


int RAMSES_CellData_Read(RAMSES_CellData_t * cell, const char * amrdata_filename, int is_hydro, int readlmax, int swap,
                         int grav_compat)
{
    int ret = READ_CELLS_NO_ERROR;
    int read_ndim, twotondim, offset, ind, i;
    int icpu, cur_ngrids, ilevel, ivarout, ivarfile;
    int more_to_read, read_this;
    void * buf = NULL;
   	FIO_File * file = NULL;

	/* Open the F90 binary file */
	file = FIO_Open(amrdata_filename, "r", 4, swap? FIO_SWAP_BYTES : FIO_DEFAULT);
	if (file == NULL)
	{
		ret = FIO_OPEN_FAIL;
		goto cleanup;
	}

	/* Read ncpu record */
	if (FIO_ReadRecord(&(cell->ncpu), 4, 1, file) != 1)
	{
		fprintf(stderr, "[DEBUG] Error while trying to read 'ncpu' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	if (!is_hydro) { /* Gravity file => no 'nvar_file' F90 record */
		/* Read nvar_file record only */
		if (FIO_ReadRecord(&(cell->nvar_file), 4, 1, file) != 1)
		{
			fprintf(stderr, "[DEBUG] Error while trying to read 'nvar' record.\n");
			ret = FIO_READ_RECORD_FAIL;
			goto cleanup;
		}
		if (grav_compat)
			cell->nvar_file++;
	} else {
		/* Read nvars and ndim records */
		if (FIO_ReadRecord(&(cell->nvar_file), 4, 1, file) != 1 || FIO_ReadRecord(&read_ndim, 4, 1, file) != 1)
		{
			fprintf(stderr, "[DEBUG] Error while trying to read 'nvar' or 'ndim' record.\n");
			ret = FIO_READ_RECORD_FAIL;
			goto cleanup;
		}
	}

	if (is_hydro && cell->ndim != read_ndim)
	{
		ret = READ_CELLS_NDIM_INVALID;
		goto cleanup;
	}

	twotondim = 1 << (cell->ndim);

	/* Check the variables index value validity */
	for(i=0; i<cell->nvar_data; i++)
	{
		if (cell->ivar_arr[i] >= cell->nvar_file)
		{
			ret = READ_CELLS_UNAVAIL_IDATA_ERROR + cell->ivar_arr[i];
			goto cleanup;
		}
	}

	/* Read levelmax and nboundary records */
	if (FIO_ReadRecord(&(cell->levelmax), 4, 1, file) != 1 || FIO_ReadRecord(&(cell->nboundary), 4, 1, file) != 1)
	{
		fprintf(stderr, "[DEBUG] Error while trying to read 'levelmax' or 'nboundary' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	 /* Skip 'gamma' record only for files other than gravity files */
	if (is_hydro)
	{
		if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL)
		{
			fprintf(stderr, "[DEBUG] Error while trying to skip 'gamma' record.\n");
			ret = FIO_SKIP_RECORD_FAIL;
			goto cleanup;
		}
	}

	/* Allocate buffer to read variable values */
	buf = malloc(sizeof(double)*cell->ngrids);
	if (buf == NULL)
	{
	    ret = READ_CELLS_MEMORY_ERROR;
	    goto cleanup;
	}

	/* Fill cell-based 3D data array : cell->data[ivar, igrid, icell] */
	offset = 0;
	for (ilevel=0; ilevel<readlmax; ilevel++) {


		for (icpu=0; icpu<(cell->ncpu)+(cell->nboundary); icpu++) {

			/* Skip ilevel record */
			if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL)
			{
				fprintf(stderr, "[DEBUG] Error while trying to skip 'ilevel' record.\n");
				ret = FIO_SKIP_RECORD_FAIL;
				goto cleanup;
			}

			//printf("ilevel = %d/%d ; icpu = %d/%d\n", ilevel, readlmax-1, icpu, cell->ncpu - 1);

			/* Read ncache record */
			if (FIO_ReadRecord(&cur_ngrids, 4, 1, file) != 1)
			{
				fprintf(stderr, "[DEBUG] Error while trying to read 'ncache' record.\n");
				ret = FIO_READ_RECORD_FAIL;
				goto cleanup;
			}

			if(cur_ngrids <= 0) continue;

			if (offset+cur_ngrids > cell->ngrids)
			{
				ret = READ_CELLS_READ_NGRID_OVERFLOW;
				goto cleanup;
			}

			for (ind=0; ind<twotondim; ind++) {
				ivarout = 0;
				for (ivarfile=0; ivarfile<cell->nvar_file; ivarfile++) {

					/* Skip variable values : nothing to read anymore */
					more_to_read = (ivarout < cell->nvar_data);
					if (!more_to_read) {
						if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL)
						{
							fprintf(stderr, "[DEBUG] Error while trying to skip amr field variable #%d record.\n", ivarfile);
							ret = FIO_SKIP_RECORD_FAIL;
							goto cleanup;
						}
						continue;
					}

					/* Skip unrequired variable values */
					read_this = (ivarfile == cell->ivar_arr[ivarout]);
					if (!read_this) {
						if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL)
						{
							fprintf(stderr, "[DEBUG] Error while trying to skip amr field variable #%d record.\n", ivarfile);
							ret = FIO_SKIP_RECORD_FAIL;
							goto cleanup;
						}
						continue;
					}

					/* Read variable values from file */
					if (FIO_ReadRecord(buf, 8, cur_ngrids, file) != cur_ngrids)
					{
						fprintf(stderr, "[DEBUG] Error while trying to read amr field variable #%d record.\n", ivarfile);
						ret = FIO_READ_RECORD_FAIL;
						goto cleanup;
					}

					/* Dispatch the values in the proper data array */
					memcpy((void*)&(CELEM3D(cell->data, ivarout, ind, offset, cell->nvar_data, twotondim, cell->ngrids)),
					       buf, cur_ngrids*sizeof(double));

					ivarout += 1;
				}
			}

			offset += cur_ngrids;
		}
	}

	if (offset != cell->ngrids)
		ret = READ_CELLS_READ_NGRID_MISS;

    cleanup:

    FREE_IF_NOTNULL(buf);

	if (file != NULL)
		FIO_Close(file);

	return ret;
}



RAMSES_CellHydroHeader_t * RAMSES_CellHydro_Header_New(void)
{
	RAMSES_CellHydroHeader_t * hdr;
	hdr = malloc(sizeof(RAMSES_CellHydroHeader_t));
	if (hdr != NULL)
	{
		hdr->nvar_file = 0;
		hdr->ncpu = 0;
		hdr->ndim = 0;
		hdr->levelmax = 0;
		hdr->gamma = 0.0;
	}
	return hdr;
}


void RAMSES_CellHydro_Header_Free(RAMSES_CellHydroHeader_t * hdr)
{
	FREE_IF_NOTNULL(hdr);
}


int RAMSES_CellHydro_Header_Read(RAMSES_CellHydroHeader_t * header,	const char * hydrodata_filename, int swap)
{
	int ret = READ_CELLS_NO_ERROR;
	FIO_File * file = NULL;

	/* Open the F90 binary file */
	file = FIO_Open(hydrodata_filename, "r", 4, swap? FIO_SWAP_BYTES : FIO_DEFAULT);
	if (file == NULL)
	{
		ret = FIO_OPEN_FAIL;
		goto cleanup;
	}

	/* Read ncpu record */
	if (FIO_ReadRecord(&(header->ncpu), 4, 1, file) != 1)
	{
		fprintf(stderr, "[DEBUG] Error while trying to read 'ncpu' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Read nvars and ndim records */
	if (FIO_ReadRecord(&(header->nvar_file), 4, 1, file) != 1 || FIO_ReadRecord(&header->ndim, 4, 1, file) != 1)
	{
		fprintf(stderr, "[DEBUG] Error while trying to read 'nvar' or 'ndim' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Read levelmax and nboundary records */
	if (FIO_ReadRecord(&(header->levelmax), 4, 1, file) != 1)
	{
		fprintf(stderr, "[DEBUG] Error while trying to skip 'levelmax' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Skip nboundary record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL)
	{
		fprintf(stderr, "[DEBUG] Error while trying to skip 'nboundary' record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	/* Read 'gamma' record */
	if (FIO_ReadRecord(&(header->gamma), 8, 1, file) != 1)
	{
		fprintf(stderr, "[DEBUG] Error while trying to read 'gamma' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	cleanup:

	if (file != NULL)
		FIO_Close(file);

	return ret;
}

