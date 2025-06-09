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
#include "read_amr.h"


RAMSES_AmrHeader_t * RAMSES_AmrHeader_New() {
	RAMSES_AmrHeader_t * hdr = malloc(sizeof(RAMSES_AmrHeader_t));
	if (hdr != NULL)
	{
		hdr->out_times = NULL;
		hdr->out_aexps = NULL;
	}
	return hdr;
}

void RAMSES_AmrHeader_Free(RAMSES_AmrHeader_t * header) {
    if (header != NULL) {
	    FREE_IF_NOTNULL(header->out_times);
	    FREE_IF_NOTNULL(header->out_aexps);
	    free(header);
	}
}

RAMSES_AmrStruct_t * RAMSES_AmrStruct_New() {
	RAMSES_AmrStruct_t * amr = malloc(sizeof(RAMSES_AmrStruct_t));
	if (amr != NULL)
	{
		amr->ngridlevel = NULL;
		amr->ngridbound = NULL;
		amr->grid_indices = NULL;
		amr->coarse_son_indices = NULL;
		amr->grid_centers = NULL;
		amr->son_indices = NULL;
	}
	return amr;
}

void RAMSES_AmrStruct_Free(RAMSES_AmrStruct_t * amr) {
    if (amr != NULL) {
	    FREE_IF_NOTNULL(amr->ngridlevel);
	    FREE_IF_NOTNULL(amr->ngridbound);
	    FREE_IF_NOTNULL(amr->grid_indices);
	    FREE_IF_NOTNULL(amr->coarse_son_indices);
	    FREE_IF_NOTNULL(amr->grid_centers);
	    FREE_IF_NOTNULL(amr->son_indices);
	    free(amr);
	}
}

int RAMSES_AmrHeader_Read(RAMSES_AmrHeader_t * header, FIO_File * file) {
    int ret = READ_AMR_NO_ERROR;
    int i;
    int output_info[3];
    double aexp_state[5];

	/* Read ncpu record */
	if (FIO_ReadRecord(&(header->ncpu), 4, 1, file) != 1) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'ncpu' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

    /* Read ndim record */
	if (FIO_ReadRecord(&(header->ndim), 4, 1, file) != 1) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'ndim' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Read coarse grid shape record */
	if (FIO_ReadRecord(&(header->coarse_shape), 4, 3, file) != 3) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'nx,ny,nz' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	header->ncoarse = 1;
	for (i=0; i<3; i++)
		header->ncoarse *= header->coarse_shape[i];

    /* Read levelmax record */
	if (FIO_ReadRecord(&(header->levelmax), 4, 1, file) != 1) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'levelmax' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

    /* Read ngridmax record */
	if (FIO_ReadRecord(&(header->ngridmax), 4, 1, file) != 1) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'levelmax' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

    /* Read nboundary record */
	if (FIO_ReadRecord(&(header->nboundary), 4, 1, file) != 1) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'nboundary' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

    /* Read ngrids_lmax record */
	if (FIO_ReadRecord(&(header->ngrids_lmax), 4, 1, file) != 1) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'ngrids_lmax' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

    /* Read boxlen record */
	if (FIO_ReadRecord(&(header->boxlen), 8, 1, file) != 1) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'boxlen' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

    /* Read output_info record */
	if (FIO_ReadRecord(output_info, 4, 3, file) != 3) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'output_info' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	header->noutputs = output_info[0];
	header->iout = output_info[1];

	header->out_times = malloc(sizeof(double)*header->noutputs);
	if (header->out_times == NULL) {
	    ret = READ_AMR_MEMORY_ERROR;
	    goto cleanup;
	}

    /* Read out_times record */
	if (FIO_ReadRecord(header->out_times, 8, header->noutputs, file) != header->noutputs) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'out_times' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	header->out_aexps = malloc(sizeof(double)*header->noutputs);
	if (header->out_aexps == NULL) {
	    ret = READ_AMR_MEMORY_ERROR;
	    goto cleanup;
	}

	/* Read out_aexps record */
	if (FIO_ReadRecord(header->out_aexps, 8, header->noutputs, file) != header->noutputs) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'out_aexps' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Read time record */
	if (FIO_ReadRecord(&(header->time), 8, 1, file) != 1) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'time' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Skip dt_old record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip 'dt_old' record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	/* Skip dt_new record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip 'dt_new' record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	/* Skip (nstep, nstep_coarse) record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip (nstep, nstep_coarse) record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	/* Skip (const, mass_tot_0, rho_tot) record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip (const, mass_tot_0, rho_tot) record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	/* Skip cosmo parameters record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip cosmo parameters record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	/* Read aexp_state record */
	if (FIO_ReadRecord(aexp_state, 8, 5, file) != 5) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'aexp_state' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}
	header->aexp = aexp_state[0];

	/* Skip 'mass_sph' record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip 'mass_sph' record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	return ret;

	/* Cleanup */
	cleanup:

	return ret;
}

int RAMSES_ReadAmr(RAMSES_AmrHeader_t * header, RAMSES_AmrStruct_t * amr,
		           const char * amr_filename, int readlmax, int swap) {
	FIO_File * file = NULL;
	int i, j, icell, idim, igrid, ilevel, icpu, cur_ngrids, ngrids;
	int ind, offset, son, ret = READ_AMR_NO_ERROR;
//	char ordering[128];
	int ndim, twotondim, ncpu, nboundary, levelmax, *tmp = NULL, xbound[3];
    void * gbuf = NULL;

	/* Open the file */
	file = FIO_Open(amr_filename, "r", 4, swap? FIO_SWAP_BYTES : FIO_DEFAULT);
	if (file == NULL)
	{
		ret = FIO_OPEN_FAIL;
		goto cleanup;
	}

	/* Read the header and fill some of the AMR structure */
	ret = RAMSES_AmrHeader_Read(header, file);
	if (ret != READ_AMR_NO_ERROR) goto cleanup;

	amr->ndim = header->ndim;
	amr->twotondim = 1<<amr->ndim;
	amr->ncpu = header->ncpu;
	amr->nboundary = header->nboundary;
	amr->levelmax = header->levelmax;
	amr->ncoarse = header->ncoarse;
	amr->readlmax = MIN(readlmax, amr->levelmax);

	ndim = amr->ndim;
	twotondim = amr->twotondim;
	ncpu = amr->ncpu;
	nboundary = amr->nboundary;
	levelmax = amr->levelmax;

	/* Skip Grid linked lists */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip 'headl' record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip 'taill' record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	/* Grids by level and CPU */
	tmp = malloc(sizeof(int)*ncpu*levelmax);
	if (tmp == NULL) {
	    ret = READ_AMR_MEMORY_ERROR;
	    goto cleanup;
	}
	if (FIO_ReadRecord(tmp, 4, ncpu*levelmax, file) != ncpu*levelmax) {
		fprintf(stderr, "[DEBUG] Error while trying to read 'ngridlevel' record.\n");
		ret = FIO_READ_RECORD_FAIL;
		goto cleanup;
	}

	/* Remap to C order */
	amr->ngridlevel = malloc(sizeof(int)*ncpu*levelmax);
	if (amr->ngridlevel == NULL) {
	    ret = READ_AMR_MEMORY_ERROR;
	    goto cleanup;
	}
	for (i=0; i<ncpu; i++)
		for (j=0; j<levelmax; j++)
			CELEM2D(amr->ngridlevel, i, j, ncpu, levelmax) = FELEM2D(tmp, i, j, ncpu, levelmax);
	FREE_IF_NOTNULL(tmp);
    tmp = NULL;

    /* Skip 'numbtot' record */
	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
		fprintf(stderr, "[DEBUG] Error while trying to skip 'numbtot' record.\n");
		ret = FIO_SKIP_RECORD_FAIL;
		goto cleanup;
	}

	/* Boundary info */
	if (nboundary > 0) {
		/* Skip boundary grid head/tail lists */
	    if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
    		fprintf(stderr, "[DEBUG] Error while trying to skip 'headb' record.\n");
    		ret = FIO_SKIP_RECORD_FAIL;
    		goto cleanup;
    	}
    	if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
    		fprintf(stderr, "[DEBUG] Error while trying to skip 'tailb' record.\n");
    		ret = FIO_SKIP_RECORD_FAIL;
    		goto cleanup;
    	}

		tmp = malloc(sizeof(int)*nboundary*levelmax);
	    if (tmp == NULL) {
    	    ret = READ_AMR_MEMORY_ERROR;
    	    goto cleanup;
    	}
    	if (FIO_ReadRecord(tmp, 4, nboundary*levelmax, file) != nboundary*levelmax) {
	    	fprintf(stderr, "[DEBUG] Error while trying to read 'numbb' record.\n");
    		ret = FIO_READ_RECORD_FAIL;
    		goto cleanup;
    	}

		/* Remap to C order */
		amr->ngridbound = malloc(sizeof(int)*nboundary*levelmax);
		if (amr->ngridbound == NULL) {
    	    ret = READ_AMR_MEMORY_ERROR;
    	    goto cleanup;
    	}
		for (i=0; i<nboundary; i++)
			for (j=0; j<levelmax; j++)
				CELEM2D(amr->ngridbound, i, j, nboundary, levelmax) = FELEM2D(tmp, i, j, nboundary, levelmax);
		FREE_IF_NOTNULL(tmp);

		for(i=0; i<3; i++)
			xbound[i] = header->coarse_shape[i]/2;
	} else {
		amr->ngridbound = NULL;
		for(i=0; i<3; i++)
			xbound[i] = 0;
	}

    /* Skip free memory record */
    if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
        fprintf(stderr, "[DEBUG] Error while trying to skip 'free memory' record.\n");
    	ret = FIO_SKIP_RECORD_FAIL;
    	goto cleanup;
    }

    /* Ordering type */
    /* Skip ordering type record */
    if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
        fprintf(stderr, "[DEBUG] Error while trying to skip 'bisec_wall' record.\n");
    	ret = FIO_SKIP_RECORD_FAIL;
    	goto cleanup;
    }
//    /* Read ordering type */
//    if (FIO_ReadRecord(ordering, 1, 128, file) != 128) {
//	    fprintf(stderr, "[DEBUG] Error while trying to read 'ordering' record.\n");
//    	ret = FIO_READ_RECORD_FAIL;
//    	goto cleanup;
//    }

//    if (strcmp(ordering, "bisection") == 0) {
//        /* Skip bisec_wall record */
//        if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
//            fprintf(stderr, "[DEBUG] Error while trying to skip 'bisec_wall' record.\n");
//    	    ret = FIO_SKIP_RECORD_FAIL;
//    	    goto cleanup;
//        }
//
//        /* Skip bisec_next record */
//        if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
//            fprintf(stderr, "[DEBUG] Error while trying to skip 'bisec_next' record.\n");
//    	    ret = FIO_SKIP_RECORD_FAIL;
//    	    goto cleanup;
//        }
//
//        /* Skip bisec_indx record */
//        if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
//            fprintf(stderr, "[DEBUG] Error while trying to skip 'bisec_indx' record.\n");
//    	    ret = FIO_SKIP_RECORD_FAIL;
//    	    goto cleanup;
//        }
//
//        /* Skip bisec_cpubox_min record */
//        if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
//            fprintf(stderr, "[DEBUG] Error while trying to skip 'bisec_cpubox_min' record.\n");
//    	    ret = FIO_SKIP_RECORD_FAIL;
//    	    goto cleanup;
//        }
//
//        /* Skip bisec_cpubox_max record */
//        if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
//            fprintf(stderr, "[DEBUG] Error while trying to skip 'bisec_cpubox_max' record.\n");
//    	    ret = FIO_SKIP_RECORD_FAIL;
//    	    goto cleanup;
//        }
//
//    } else {
        /* Hilbert ordering */
        /* Skip Hilbert bound keys record */
        if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
            fprintf(stderr, "[DEBUG] Error while trying to skip 'bound_keys' record.\n");
    	    ret = FIO_SKIP_RECORD_FAIL;
    	    goto cleanup;
        }
//    }

	/* Cells at coarse level (level 0) */
	amr->coarse_son_indices = malloc(sizeof(int)*amr->ncoarse);
	if (amr->coarse_son_indices == NULL) {
	    ret = READ_AMR_MEMORY_ERROR;
    	goto cleanup;
    }
    if (FIO_ReadRecord(amr->coarse_son_indices, 4, amr->ncoarse, file) != amr->ncoarse) {
	    fprintf(stderr, "[DEBUG] Error while trying to read 'coarse_son_indices' record.\n");
    	ret = FIO_READ_RECORD_FAIL;
    	goto cleanup;
    }
	for (icell=0; icell<amr->ncoarse; icell++) amr->coarse_son_indices[icell] -= 1;

    /* Skip Coarse flags record */
    if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
        fprintf(stderr, "[DEBUG] Error while trying to skip 'flag' record.\n");
        ret = FIO_SKIP_RECORD_FAIL;
        goto cleanup;
    }
    /* Skip Coarse cpumap record */
    if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
        fprintf(stderr, "[DEBUG] Error while trying to skip 'cpumap' record.\n");
        ret = FIO_SKIP_RECORD_FAIL;
        goto cleanup;
    }

	/* Precompute total number of fine grids */
	ngrids = 0;
	for (ilevel=0; ilevel<amr->readlmax; ilevel++) {
		for (icpu=0; icpu<ncpu+nboundary; icpu++)
			if(icpu<ncpu)
				ngrids += CELEM2D(amr->ngridlevel, icpu, ilevel, ncpu, levelmax);
			else
				ngrids += CELEM2D(amr->ngridbound, icpu-ncpu, ilevel, nboundary, levelmax);
	}
	amr->ngrids = ngrids;

	/* Temporary buffer for holding grid data records. Size is at most 8*ngrids,
	 * if data is written in double precision. */
	gbuf = malloc(8*ngrids);
	amr->grid_indices = malloc(sizeof(int)*ngrids);
	amr->grid_centers = malloc(sizeof(double)*ngrids*ndim);
	amr->son_indices  = malloc(sizeof(int)*ngrids*twotondim);

	if (gbuf == NULL || amr->grid_indices == NULL || amr->grid_centers == NULL || amr->son_indices == NULL) {
	    ret = READ_AMR_MEMORY_ERROR;
    	goto cleanup;
	}

	/* Loop over fine levels */
	offset = 0;
	for (ilevel=0; ilevel<amr->readlmax; ilevel++) {
		for (icpu=0; icpu<ncpu+nboundary; icpu++) {
			if (icpu<ncpu)
				cur_ngrids = CELEM2D(amr->ngridlevel, icpu, ilevel, ncpu, levelmax);
			else
				cur_ngrids = CELEM2D(amr->ngridbound, icpu-ncpu, ilevel, nboundary, levelmax);

			if(cur_ngrids <= 0) continue;

			/* Grid index */
	        if (FIO_ReadRecord(gbuf, 4, cur_ngrids, file) != cur_ngrids) {
	            fprintf(stderr, "[DEBUG] Error while trying to read 'ind_grid' record.\n");
	            ret = FIO_READ_RECORD_FAIL;
	            goto cleanup;
	        }
			for (igrid=0; igrid<cur_ngrids; igrid++)
				amr->grid_indices[igrid+offset] = ((int*)gbuf)[igrid] - 1;

			/* Linked lists: prev and next */
            if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
                fprintf(stderr, "[DEBUG] Error while trying to skip 'next(ind_grid)' record.\n");
                ret = FIO_SKIP_RECORD_FAIL;
                goto cleanup;
            }
            if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
                fprintf(stderr, "[DEBUG] Error while trying to skip 'prev(ind_grid)' record.\n");
                ret = FIO_SKIP_RECORD_FAIL;
                goto cleanup;
            }
			
			/* Grid centers */
			for (idim=0; idim<ndim; idim++) {
			    if (FIO_ReadRecord(gbuf, 8, cur_ngrids, file) != cur_ngrids) {
	                fprintf(stderr, "[DEBUG] Error while trying to read 'xg(ind_grid)' record.\n");
    	            ret = FIO_READ_RECORD_FAIL;
    	            goto cleanup;
    	        }
				for (igrid=0; igrid<cur_ngrids; igrid++)
					CELEM2D(amr->grid_centers, igrid+offset, idim, ngrids, ndim) = ((double*)gbuf)[igrid] - xbound[idim];
			}

			/* Father index */
			if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
                fprintf(stderr, "[DEBUG] Error while trying to skip 'father(ind_grid)' record.\n");
                ret = FIO_SKIP_RECORD_FAIL;
                goto cleanup;
            }

			/* Neighbor index */
			for (ind=0; ind<2*ndim; ind++) {
    			if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
                    fprintf(stderr, "[DEBUG] Error while trying to skip 'nbor(ind_grid)' record.\n");
                    ret = FIO_SKIP_RECORD_FAIL;
                    goto cleanup;
                }
			}

			/* Son index */
			for (ind=0; ind<twotondim; ind++) {
			    if (FIO_ReadRecord(gbuf, 4, cur_ngrids, file) != cur_ngrids) {
	                fprintf(stderr, "[DEBUG] Error while trying to read son(ind_grid)' record.\n");
    	            ret = FIO_READ_RECORD_FAIL;
    	            goto cleanup;
    	        }
				for (igrid=0; igrid<cur_ngrids; igrid++) {
					son = ((int*)gbuf)[igrid];
					/* Remap son index */
					son -= 1;
					CELEM2D(amr->son_indices, igrid+offset, ind, ngrids, twotondim) = son;
				}
			}
				
			/* Cell CPU map */
			for(ind=0; ind<twotondim; ind++) {
    			if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
                    fprintf(stderr, "[DEBUG] Error while trying to skip 'cpu_map(ind_grid)' record.\n");
                    ret = FIO_SKIP_RECORD_FAIL;
                    goto cleanup;
                }
			}

			/* Cell refinement map */
			for(ind=0; ind<twotondim; ind++) {
    			if (FIO_SkipRecord(file) == FIO_SKIP_RECORD_FAIL) {
                    fprintf(stderr, "[DEBUG] Error while trying to skip 'flag(ind_grid)' record.\n");
                    ret = FIO_SKIP_RECORD_FAIL;
                    goto cleanup;
                }
			}

			offset += cur_ngrids;
		}
	}

    if (offset != ngrids) {
        fprintf(stderr, "[DEBUG] Number of AMR grids read mismatch.\n");
        ret = RAED_AMR_NGRID_MISMATCH_ERROR;
        goto cleanup;
    }

	/* Cleanup */
	FREE_IF_NOTNULL(gbuf);
	FREE_IF_NOTNULL(tmp);

	FIO_Close(file);

	return READ_AMR_NO_ERROR;

	/* Cleanup */
	cleanup:

	FREE_IF_NOTNULL(gbuf);
	FREE_IF_NOTNULL(tmp);

	if (file != NULL) FIO_Close(file);

	return ret;
}
