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
#include <errno.h>
#include "byteswap.h"
#include "fio.h"

FIO_File * FIO_Open(const char * filename, const char * mode, size_t markersize, int flags) {
	FIO_File * f = NULL;

	if (markersize != sizeof(int))
	{
		fprintf(stderr, "[FIO_Open] 'markersize' parameter does not have a valid size_t value");
		return NULL;
	}

	f = malloc(sizeof(FIO_File));
	if (f != NULL)
	{
		f->file = fopen(filename, mode);
		if (f->file == NULL) {
			fprintf(stderr, "[FIO_Open] Failed to open file.\n");
			free(f);
			return NULL;
		}
		f->markersize = markersize;
		f->flags = flags;
	}
	else
	    fprintf(stderr, "[FIO_Open] Could not allocate memory for a FIO_File.\n");
	return f;
}

int FIO_Close(FIO_File * fiofile) {
	int ierr;

	if (fiofile == NULL)
	{
		fprintf(stderr, "[FIO_Close] File is NULL.\n");
		return FIO_CLOSE_FAIL;
	}

	ierr = fclose(fiofile->file);
	if (ierr) {
		fprintf(stderr, "[FIO_Close] Error while closing file.\n");
		free(fiofile);
		fiofile = NULL;
		return FIO_CLOSE_FAIL;
	}

	free(fiofile);
	fiofile = NULL;

	return ierr;
}

/* Low-level read function that handles byteswapping */
size_t _FIO_FRead(FIO_File * fiofile, void * buf, size_t itemsize, size_t count) {
	if (fiofile == NULL || fiofile->file == NULL)
	{
		fprintf(stderr, "[_FIO_FRead] File is NULL.\n");
		return (size_t)0;
	}

    if (itemsize == 0 || count == 0)
    {
        return (size_t)0;
    }

	size_t nread = fread(buf, itemsize, count, fiofile->file);

	if(fiofile->flags & FIO_SWAP_BYTES) {
		size_t i;
		switch(itemsize) {
			case 4:
				for(i=0; i<nread; i++) {
					unsigned int *ptr = (unsigned int*)buf + i;
					*ptr = __bswap_32(*ptr);
				}
				break;

			case 8:
				for(i=0; i<nread; i++) {
					unsigned long *ptr = (unsigned long*)buf + i;
					*ptr = __bswap_64(*ptr);
				}
				break;

			default:
				fprintf(stderr, "[_FIO_FRead] Invalid item size\n");
				return (size_t)0;
		}
	}

	return nread;
}

int _FIO_ReadMarker(FIO_File * fiofile) {
	size_t nread;
	//int32_t nread;
	size_t marker;
	//int32_t marker;  // Ensure 4-byte marker
	//size_t nread = _FIO_FRead(fiofile, &marker, sizeof(int32_t), 1);
	
	size_t *ptr = &marker;
	

	nread = _FIO_FRead(fiofile, ptr, fiofile->markersize, 1);

	if (nread == 0)
	{
		fprintf(stderr, "[_FIO_ReadMarker] No data was read.\n");
		return 0;
	}

	switch (fiofile->markersize) {
		case 4:
			return (int)(*((unsigned int*)ptr));
			break;
		case 8:
			return (int)(*((unsigned long*)ptr));
			break;
		default:
			fprintf(stderr, "[_FIO_ReadMarker] Invalid marker size\n");
			return -1;
			break;
	}
}

int FIO_SkipRecord(FIO_File * fiofile) {
	int m1, m2;

	m1 = _FIO_ReadMarker(fiofile);

	if (m1 < 0)
		return FIO_SKIP_RECORD_FAIL;

	if (m1 == 0) // No data read => return 0
                return 0;

	int ierr = fseek(fiofile->file, m1, SEEK_CUR);
	if (ierr) {
		fprintf(stderr, "[FIO_SkipRecord] Seek error\n");
		return FIO_SKIP_RECORD_FAIL;
	}
	m2 = _FIO_ReadMarker(fiofile);
	if (m1 != m2) {
		fprintf(stderr, "[FIO_SkipRecord] Marker mismatch\n");
		return FIO_SKIP_RECORD_FAIL;
	}

	/* No error */
	return 0;
}


int FIO_ReadRecord(void * buffer, size_t itemsize, size_t maxcount, FIO_File * fiofile) {
	size_t m1, m2;
	m1 = _FIO_ReadMarker(fiofile);

	if (m1 < 0)
	{
		fprintf(stderr, "[FIO_ReadRecord] Corrupted F90 binary record.\n");
		return -1;
	}

        if (m1 == 0)
                return 0;

	if (m1 != itemsize*maxcount) {
		fprintf(stderr, "[FIO_ReadRecord] Read overflow : trying to read %d item(s) of %d bytes but found F90 record marker "
		"of %d bytes.\n", (int)maxcount, (int)itemsize, (int)m1);
		return -1;
	}

	size_t nread = _FIO_FRead(fiofile, buffer, itemsize, maxcount);
	if(nread != maxcount) {
		fprintf(stderr, "[FIO_ReadRecord] Wrong number of items read.\n");
		return -1;
	}

	m2 = _FIO_ReadMarker(fiofile);
	if (m1 != m2) {
		fprintf(stderr, "[FIO_ReadRecord] Marker mismatch (%d != %d).\n", (int)m1, (int)m2);
		return -1;
	}

	return (int)nread;
}
