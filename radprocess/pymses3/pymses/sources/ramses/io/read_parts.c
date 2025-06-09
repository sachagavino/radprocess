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
#include <stdint.h>
#include "fio.h"
#include "utils.h"
#include "read_parts.h"


RAMSES_PartData_t * RAMSES_PartData_New()
{
	RAMSES_PartData_t * part = malloc(sizeof(RAMSES_PartData_t));
	if (part != NULL)
	{
		part->pos = NULL;
		part->vel = NULL;
		part->mass = NULL;
		part->id_32 = NULL;
		part->id_64 = NULL;
		part->level = NULL;
		part->epoch = NULL;
		part->metal = NULL;
	}
	return part;
}

void RAMSES_PartData_Free(RAMSES_PartData_t * part)
{
	FREE_IF_NOTNULL(part->pos);
	FREE_IF_NOTNULL(part->vel);
	FREE_IF_NOTNULL(part->mass);
	FREE_IF_NOTNULL(part->id_32);
	FREE_IF_NOTNULL(part->id_64);
	FREE_IF_NOTNULL(part->level);
	FREE_IF_NOTNULL(part->epoch);
	FREE_IF_NOTNULL(part->metal);
	free(part);
	part = NULL;
}


int RAMSES_PartData_Read(RAMSES_PartData_t * part, double boxlen, FIO_File * file, int long_int)
{
	/* 'long_int' parameter is 1 <=> RAMSES was compiled with the -DLONGINT flag */
    int ret = READ_PARTS_NO_ERROR;
	int read_length;
	int npart, ndim, idim, ipart;
	double * dtmp = NULL;

	/* File header */
	// Read ncpu, ndim and npart records
	if (FIO_ReadRecord(&(part->ncpu), 4, 1, file) != 1 || FIO_ReadRecord(&(part->ndim), 4, 1, file) != 1 ||
		FIO_ReadRecord(&(part->npart), 4, 1, file) != 1)
	{
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Skip localseed record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL)
	{
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	// Read nstar record
	if (FIO_ReadRecord(&(part->nstar), 4, 1, file) != 1)
	{
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Skip mstar and mstar_lost records */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL || FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL)
	{
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	// Read nsink record
	if (FIO_ReadRecord(&(part->nsink), 4, 1, file) != 1)
	{
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	npart = part->npart;
	ndim = part->ndim;
	dtmp = malloc(sizeof(double)*npart);
	if (dtmp == NULL)
	{
		ret = READ_PARTS_MEMORY_ERROR;
		goto cleanup;
	}

	/* Position */
	part->pos = malloc(sizeof(double)*npart*ndim);
	if (part->pos == NULL)
	{
		ret = READ_PARTS_MEMORY_ERROR;
		goto cleanup;
	}

	for(idim=0; idim<ndim; idim++) {
		if (FIO_ReadRecord(dtmp, 8, npart, file) != npart)
		{
			ret = FIO_READ_RECORD_FAIL;
			goto cleanup;
		}

		for(ipart=0; ipart<npart; ipart++)
			CELEM2D(part->pos, ipart, idim, npart, ndim) = dtmp[ipart]/boxlen;
	}

	/* Velocity */
	part->vel = malloc(sizeof(double)*npart*ndim);
	if (part->vel == NULL)
	{
		ret = READ_PARTS_MEMORY_ERROR;
		goto cleanup;
	}
	for(idim=0; idim<ndim; idim++) {
		if (FIO_ReadRecord(dtmp, 8, npart, file) != npart)
		{
			ret = FIO_READ_RECORD_FAIL;
			goto cleanup;
		}

		for(ipart=0; ipart<npart; ipart++)
			CELEM2D(part->vel, ipart, idim, npart, ndim) = dtmp[ipart];
	}

	/* Mass */
	part->mass = malloc(sizeof(double)*npart);
	if (part->mass == NULL)
	{
		ret = READ_PARTS_MEMORY_ERROR;
		goto cleanup;
	}

	if (FIO_ReadRecord(part->mass, 8, npart, file) != npart)
	{
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;

	}

	/* Id */
	if (long_int)
		part->id_64 = malloc(sizeof(int64_t)*npart);
	else
		part->id_32 = malloc(sizeof(int)*npart);

	if (part->id_64 == NULL && part->id_32 == NULL)
	{
		ret = READ_PARTS_MEMORY_ERROR;
		goto cleanup;
	}

	if (long_int)
	{
		if (FIO_ReadRecord(part->id_64, 8, npart, file) != npart)
		{
			ret = FIO_READ_RECORD_FAIL;
			goto cleanup;
		}
	}
	else
	{
		if (FIO_ReadRecord(part->id_32, 4, npart, file) != npart)
		{
			ret = FIO_READ_RECORD_FAIL;
			goto cleanup;
		}
	}


	/* Level */
	part->level = malloc(sizeof(int)*npart);
	if (part->level == NULL)
	{
		ret = READ_PARTS_MEMORY_ERROR;
		goto cleanup;
	}

	if (FIO_ReadRecord(part->level, 4, npart, file) != npart)
	{
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

//	/* Skip family and tag records */
//	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL || FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL)
//	{
//		ret = FIO_SKIP_RECORD_FAIL;
//		goto cleanup;
//	}

	/* Epoch & Metallicity */
	if(part->nstar > 0 || part->nsink > 0) {
		part->epoch = malloc(sizeof(double)*npart);
		if (part->epoch == NULL)
		{
			ret = READ_PARTS_MEMORY_ERROR;
			goto cleanup;
		}

		if (FIO_ReadRecord(part->epoch, 8, npart, file) != npart)
		{
			ret = FIO_READ_RECORD_FAIL;
			goto cleanup;
		}

		part->metal = malloc(sizeof(double)*npart);
		if (part->metal == NULL)
		{
			ret = READ_PARTS_MEMORY_ERROR;
			goto cleanup;
		}

		read_length = FIO_ReadRecord(part->metal, 8, npart, file);
		if (read_length != npart) {
		        fprintf(stderr, "No metallicity data available.\n");
			free(part->metal);
			part->metal = NULL;
		}
	}
	
	cleanup:

	/* Cleanup */
	FREE_IF_NOTNULL(dtmp);

	return ret;
}
