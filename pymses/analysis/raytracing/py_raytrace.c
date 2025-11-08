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
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "numpy/arrayobject.h"
#include <strings.h>

#define FREE_IF_NOTNULL(A)   if(A != NULL) free(A);

#define RAYTRACE_NO_ERROR      0
#define RAYTRACE_MEMORY_ERROR  1

/* Static pointer to a raytracing error object */
static PyObject *_RayTracingError;

/* Module docstring */
static char module_docstring[] = "This module provides RAMSES octree ray-casting functionalities.";

/* ------------------------------------------------------------------------------------------------------------------ *
 * --------- raytrace_amr() routine to compute 2D ray-traced maps from Ramses hydrodynamical 3D AMR data ------------ *
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject *
raytrace_amr(PyObject *self, PyObject *args)
{
	int ierr = RAYTRACE_NO_ERROR, jray, jray2, cell_index;
	int rlev, nblocks, nrays, nscal, max_op, use_dx;
	PyArrayObject *maps_n, *ray_l_cumul_n, *ray_orig_n, *ray_vects_n, *ray_length_max_n;
	PyArrayObject *cell_centers_n, *sons_n, *scal_data_n, *igrids_n, *icells_n, *active_mask_n, *grid_levels_n, *params_n;
	double *maps, *ray_l_cumul, *ray_orig, *ray_vects, *ray_depth_max, *cell_centers, *scal_data;
	int *sons, *igrids, *grid_levels, *active_mask, *params;
	char *icells;

	int twotondim, twotondimMndim, ndim, twotondimMnscal;
	int igrid_root, level, iblock, idim, idim2, iray, head;
	int igrid, icell, igrid_son, ilevel, iscal;
	int ind_son, intersect, sgn_face, block_already_done;
	int ind_son1, ind_son2, ind_entry, ind_exit, add, iblock2;
	char icell_root;
	double alpha, alpha_save, dx, cube_depth, r2, l, d, v, vl;
	double Rsphere, Rsphere2;
	double sqrt3 = sqrt(3.0);
	double dm2, depth_max;

	int* igrid_pile = NULL;
	int* icell_pile = NULL;
	int* level_pile = NULL;
	int* dim_constr = NULL;
	int* dim_side = NULL;
	double* block_center = NULL;

	if(!PyArg_ParseTuple(args, "OOOOOOOOOOOOO",
                         &maps_n, &ray_l_cumul_n, &ray_orig_n, &ray_vects_n, &ray_length_max_n, &cell_centers_n,
	                     &sons_n, &scal_data_n, &igrids_n, &icells_n, &active_mask_n, &grid_levels_n, &params_n))
		return NULL;

	/* Fetch the numpy array data C pointers */
	maps = PyArray_DATA(maps_n);
	ray_l_cumul = PyArray_DATA(ray_l_cumul_n);
	ray_orig = PyArray_DATA(ray_orig_n);
	ray_vects = PyArray_DATA(ray_vects_n);
	ray_depth_max = PyArray_DATA(ray_length_max_n);
	cell_centers = PyArray_DATA(cell_centers_n);
	sons = PyArray_DATA(sons_n);
	scal_data = PyArray_DATA(scal_data_n);
	igrids = PyArray_DATA(igrids_n);
	icells = PyArray_DATA(icells_n);
	grid_levels = PyArray_DATA(grid_levels_n);
	active_mask = PyArray_DATA(active_mask_n);
	params = PyArray_DATA(params_n);

	/* Get the integer parameters */
	nrays = params[0];
	ndim = params[1];
	nblocks = params[2];
	rlev = params[3];
	nscal = params[4];
	twotondim = 1 << ndim;
	twotondimMndim = twotondim * ndim;
	twotondimMnscal = twotondim * nscal;
	max_op = params[5];
	use_dx = params[6];


	igrid_pile = malloc(twotondim * rlev * sizeof(int));
	icell_pile = malloc(twotondim * rlev * sizeof(int));
	level_pile = malloc(twotondim * rlev * sizeof(int));
	dim_constr = malloc(ndim * sizeof(int));
	dim_side = malloc(ndim * sizeof(int));
	block_center = malloc(ndim * sizeof(double));

	if (igrid_pile == NULL || icell_pile == NULL || level_pile == NULL || dim_constr == NULL || dim_side == NULL ||
		block_center == NULL)
	{
		ierr = RAYTRACE_MEMORY_ERROR;
		goto cleanup;
	}

	for (iblock = 0 ; iblock < nblocks ; iblock++)
	{
		igrid_root = igrids[iblock];
		icell_root = icells[iblock];
		block_already_done = 0;
		for (iblock2 = 0 ; iblock2 < iblock ; iblock2++)
		{
			if ((igrids[iblock2] == igrids[iblock]) && (icells[iblock2] == icells[iblock])){
				block_already_done = 1;
				break;
			}
		}
		if (block_already_done == 1)
			continue;
		level = grid_levels[igrid_root];
		dx = 1.0/(1<<level);
		Rsphere = sqrt3 * dx / 2.0;
		Rsphere2 = Rsphere*Rsphere;
		for (idim = 0 ; idim < ndim ; idim++){
			block_center[idim] = cell_centers[igrid_root * twotondimMndim + icell_root * ndim + idim];
		}

		/* ################# *
		 * # Ray iteration # *
		 * ################# */
		for (iray = 0 ; iray < nrays ; iray++)
		{
			depth_max = ray_depth_max[iray];
			dm2 = depth_max/2.0;
			/* Coordinate parameter of the point of the ray nearest to the block center */
			alpha = 0.0;
			for (idim = 0 ; idim < ndim ; idim++)
			{
			    jray = iray * ndim + idim;
				alpha += (block_center[idim]-ray_orig[jray])*ray_vects[jray];
			}

			/* Square distances of the ray from the block center */
			r2 = 0.0;
			for (idim = 0 ; idim < ndim ; idim++)
			{
			    jray = iray * ndim + idim;
				v = ray_orig[jray] + alpha * ray_vects[jray] - block_center[idim];
				r2 += v*v;
			}
			/* Is this ray able to intersect the block ? */
			if (r2 > Rsphere2)// # no => continue
				continue;
			/* Distance from the ray extreme points */
			d = alpha-dm2;
			if (d<0.)
				d=-d;
			d = d-dm2;
			if (d>0.) {
				r2 += d*d;
				/* Is this ray still able to intersect the block ? */
				if (r2 > Rsphere2)
					continue;
			}

			/* Pile initialisation with current block */
			head=0;
			igrid_pile[head] = igrid_root;
			icell_pile[head] = icell_root;
			level_pile[head] = level;
			while (head !=-1)
			{
				igrid = igrid_pile[head];
				icell = icell_pile[head];
				ilevel = level_pile[head];
				head -= 1;

				/* Current cell size */
				dx = 1.0/(1<<ilevel);
				/* Cube optical depth computation */
				cube_depth = 0.0;
				intersect = 0;
				alpha_save = 0.0;
				ind_son1 = 0;
				ind_son2 = 0;
				for (idim2 = 0 ; idim2 < ndim ; idim2++)
				{
				    jray2 = iray * ndim + idim2;
					if(ray_vects[jray2] == 0.0)
					{
						/* ray is parallel to the face along dimension idim. There won't be an intersection. */
						continue;
					}

					/* All ray-cell intersections found => break loop */
					if (intersect==2) break;
					for (sgn_face = -1 ; sgn_face < 2 ; sgn_face+=2){ /* for sgn_face in [-1,1]: */
						/* Coordinate parameter of the intersection of the ray with the face plane */
						alpha = ((cell_centers[igrid * twotondimMndim + icell * ndim + idim2] + sgn_face * dx / 2.0)
							 - ray_orig[jray2])/ray_vects[jray2];
						/* Is this intersection within the face ? */
						v=0.0;
						vl=0.0;
						ind_son = (twotondim-1) - (((-sgn_face+1)>>1)<<idim2);
						for (idim = 0 ; idim < ndim ; idim++){
							if (idim != idim2)
							{
							    jray = iray * ndim + idim;
								l = (ray_orig[jray] + alpha * ray_vects[jray]
								     - cell_centers[igrid * twotondimMndim + icell * ndim + idim])/dx;
								if ((l>vl) && (ray_vects[jray] == 0.0))
									vl=l;
								if (l<0.0)
								{
									l=-l;
									ind_son -= (1<<idim);
								}
								if (l>v)
									v=l;
							}
						}
						if ((v <= 0.5) && (vl < 0.5)){
							/* The intersection is in the face of the cube */
							intersect +=1;
							if (intersect == 1){ //# First intersection found
								if (alpha < 0.0)
									alpha_save = 0.0;
								else if (alpha > depth_max)
									alpha_save = depth_max;
								else
									alpha_save = alpha;
								ind_son1=ind_son;
							}
							else{// # Second intersection found : done
								if (alpha < 0.0)
									alpha = 0.0;
								else if (alpha > depth_max)
									alpha = depth_max;
								cube_depth = alpha-alpha_save;
									/* Make sure ind_son1 = entry point son index
									             ind_son2 = exit point son index */
								if (cube_depth < 0.0) {// # ind_son is the entry point son index
									cube_depth=-cube_depth;
									ind_son2 = ind_son1;
									ind_son1 = ind_son;
								}
								else//# ind_son is the exit point son index
									ind_son2=ind_son;
								break;
							}
						}
					}
				}
				/* Is the current cell actually intersected by the ray ? */
				if (cube_depth > 0.0)
				{
					/* The ray cuts the block */
					cell_index = igrid * twotondim + icell;
					igrid_son = sons[cell_index];
					if ((igrid_son < 0) || (ilevel >= rlev))
					{/* We're done searching: the subcell is the one we want */
						// We skip inactive cells
						if (active_mask[igrid] == 0)
			                continue;

						/* Add cell content to the maps */
						for (iscal = 0 ; iscal < nscal ; iscal++){
							if (use_dx)
								v = ilevel * 1.0;
							else
								v = scal_data[igrid * twotondimMnscal + icell * nscal + iscal];
							if (max_op){
								if (v > maps[iray * nscal + iscal])
									maps[iray * nscal + iscal] = v;
							}
							else
								maps[iray * nscal + iscal] += v * cube_depth;
						}

						/* Add cube depth to the computed ray length */
						ray_l_cumul[iray] += cube_depth;
					}
					else{ /* Go down the AMR tree */
						/* ######################################################################################### *
						 * #                      Add the NECESSARY son cells to the pile                          # *
						 * ######################################################################################### */
						if (ind_son1 != ind_son2)
						{
							for (idim = 0 ; idim < ndim ; idim++){
								ind_entry = (ind_son1>>idim) & 1;
								ind_exit = (ind_son2>>idim) & 1;
								if (ind_entry==ind_exit){
									dim_constr[idim] = 1;
									dim_side[idim] = ind_entry;
								}
								else
									dim_constr[idim] = 0;
							}

							/* Add the POSSIBLY intersecting son cells different from the entry/exit point son cells */
							for (ind_son = 0 ; ind_son < twotondim ; ind_son++){ // # All the son cells
								add = 1;
								/* Different from entry/exit point son cells */
								if ((ind_son == ind_son1)||(ind_son == ind_son2))
									add = 0;
								else
								{
									for (idim = 0 ; idim < ndim ; idim++){
										if ((dim_constr[idim])&&(dim_side[idim] != ((ind_son>>idim) & 1))){
											add = 0;
											break;
										}
									}
								}

								if (add)
								{
									head += 1;
									level_pile[head] = ilevel+1;
									igrid_pile[head] = igrid_son;
									icell_pile[head] = ind_son;
								}
							}

							/* If different from the son cell of the entry point, add the exit point son cell */
							head += 1;
							level_pile[head] = ilevel+1;
							igrid_pile[head] = igrid_son;
							icell_pile[head] = ind_son2;
						}

						/* Add the son cell of the entry point of the ray in the father grid */
						head+=1;
						level_pile[head] = ilevel+1;
						igrid_pile[head] = igrid_son;
						icell_pile[head] = ind_son1;
					}
				}
			}
		}
	}

	cleanup:

	FREE_IF_NOTNULL(igrid_pile);
	FREE_IF_NOTNULL(icell_pile);
	FREE_IF_NOTNULL(level_pile);
	FREE_IF_NOTNULL(dim_constr);
	FREE_IF_NOTNULL(dim_side);
	FREE_IF_NOTNULL(block_center);

	if (ierr == RAYTRACE_MEMORY_ERROR)
	{
		PyErr_NoMemory();
		//Py_RETURN_NONE;
		return NULL;
	}

	Py_RETURN_NONE;
}


/* ------------------------------------------------------------------------------------------------------------------ *
 * ----- raytrace_amr_spherical() routine to compute 2D ray-traced maps from Ramses hydrodynamical 3D AMR data ------ *
 * ------------------------------------------------------------------------------------------------------------------ */
//#define DEBUG_PRINT 1

static PyObject *
raytrace_amr_spherical(PyObject *self, PyObject *args)
{
	int ierr = RAYTRACE_NO_ERROR, jray, jray2, cell_index;
	int rlev, nblocks, nrays, nscal, max_op, use_dx;
	PyArrayObject *maps_n, *ray_l_cumul_n, *ray_orig_n, *ray_vects_n, *ray_length_max_n;
	PyArrayObject *cell_centers_n, *sons_n, *scal_data_n, *igrids_n, *icells_n, *active_mask_n, *grid_levels_n;
	PyArrayObject *dS_sr_n, *front_clip_n, *params_n;
	double *maps, *ray_l_cumul, *ray_orig, *ray_vects, *ray_depth_max, *cell_centers, *scal_data, *dS_sr, *front_clip;
	int *sons, *igrids, *grid_levels, *active_mask, *params;
	char *icells;

	int twotondim, twotondimMndim, ndim, twotondimMnscal;
	int igrid_root, level, iblock, idim, idim2, iray, head;
	int igrid, icell, igrid_son, ilevel, iscal;
	int ind_son, intersect, sgn_face, block_already_done;
	int ind_son1, ind_son2, ind_entry, ind_exit, add, iblock2;
	char icell_root;
	double alpha, alpha_min, alpha_max, alpha_mem, dist0, dist1, dx, cube_depth, r2, l, d, v, vl;
	double Rsphere, Rsphere2;
	double sqrt3 = sqrt(3.0);
	double p1, p2, p3;
	double dm2, depth_max, ray_front_clip;

	int* igrid_pile = NULL;
	int* icell_pile = NULL;
	int* level_pile = NULL;
	int* dim_constr = NULL;
	int* dim_side = NULL;

	int fair_intersection;
	double* block_center = NULL;

	if(!PyArg_ParseTuple(args, "OOOOOOOOOOOOOOO",
                         &maps_n, &ray_l_cumul_n, &ray_orig_n, &ray_vects_n, &ray_length_max_n, &dS_sr_n, &front_clip_n,
                         &cell_centers_n, &sons_n, &scal_data_n, &igrids_n, &icells_n, &active_mask_n, &grid_levels_n,
                         &params_n))
		return NULL;

	/* Fetch the numpy array data C pointers */
	maps = PyArray_DATA(maps_n);
	ray_l_cumul = PyArray_DATA(ray_l_cumul_n);
	ray_orig = PyArray_DATA(ray_orig_n);
	ray_vects = PyArray_DATA(ray_vects_n);
	ray_depth_max = PyArray_DATA(ray_length_max_n);
	dS_sr = PyArray_DATA(dS_sr_n);
	front_clip = PyArray_DATA(front_clip_n);
	cell_centers = PyArray_DATA(cell_centers_n);
	sons = PyArray_DATA(sons_n);
	scal_data = PyArray_DATA(scal_data_n);
	igrids = PyArray_DATA(igrids_n);
	icells = PyArray_DATA(icells_n);
	grid_levels = PyArray_DATA(grid_levels_n);
	active_mask = PyArray_DATA(active_mask_n);
	params = PyArray_DATA(params_n);

	/* Get the integer parameters */
	nrays = params[0];
	ndim = params[1];
	nblocks = params[2];
	rlev = 40; //params[3];
	nscal = params[4];
	twotondim = 1 << ndim;
	twotondimMndim = twotondim * ndim;
	twotondimMnscal = twotondim * nscal;
	max_op = params[5];
	use_dx = params[6];


	igrid_pile = malloc(twotondim * rlev * sizeof(int));
	icell_pile = malloc(twotondim * rlev * sizeof(int));
	level_pile = malloc(twotondim * rlev * sizeof(int));
	dim_constr = malloc(ndim * sizeof(int));
	dim_side = malloc(ndim * sizeof(int));
	block_center = malloc(ndim * sizeof(double));

	if (igrid_pile == NULL || icell_pile == NULL || level_pile == NULL || dim_constr == NULL || dim_side == NULL ||
		block_center == NULL)
	{
		ierr = RAYTRACE_MEMORY_ERROR;
		goto cleanup;
	}

	for (iblock = 0 ; iblock < nblocks ; iblock++)
	{
		igrid_root = igrids[iblock];
		icell_root = icells[iblock];
		block_already_done = 0;
		for (iblock2 = 0 ; iblock2 < iblock ; iblock2++)
		{
			if ((igrids[iblock2] == igrids[iblock]) && (icells[iblock2] == icells[iblock])){
				block_already_done = 1;
				break;
			}
		}
		if (block_already_done == 1)
			continue;
		level = grid_levels[igrid_root];

		dx = 1.0/(1<<level);
		Rsphere = sqrt3 * dx / 2.0;
		Rsphere2 = Rsphere*Rsphere;

		for (idim = 0 ; idim < ndim ; idim++)
			block_center[idim] = cell_centers[igrid_root * twotondimMndim + icell_root * ndim + idim];

#ifdef DEBUG_PRINT
		printf("------------------------------------- New block ------------------------------\n");
		for (idim = 0 ; idim < ndim ; idim++)
			printf("block_center[%d] = %g\n", idim, block_center[idim]);
		printf("block_size = %g\n", dx);
		printf("------------------------------------------------------------------------------\n");
#endif

		/* ################# *
		 * # Ray iteration # *
		 * ################# */
		for (iray = 0 ; iray < nrays ; iray++)
		{
			depth_max = ray_depth_max[iray];
			if (front_clip[iray] > 0.0 && front_clip[iray] <= depth_max)
				ray_front_clip = depth_max - front_clip[iray];
			else
			    ray_front_clip = depth_max;

			dm2 = depth_max/2.0;
			/* Coordinate parameter of the point of the ray nearest to the block center */
			alpha = 0.0;
			for (idim = 0 ; idim < ndim ; idim++)
			{
			    jray = iray * ndim + idim;
				alpha += (block_center[idim]-ray_orig[jray])*ray_vects[jray];
			}

			/* Square distances of the ray from the block center */
			r2 = 0.0;
			for (idim = 0 ; idim < ndim ; idim++)
			{
			    jray = iray * ndim + idim;
				v = ray_orig[jray] + alpha * ray_vects[jray] - block_center[idim];
				r2 += v*v;
			}
			/* Is this ray able to intersect the block ? */
			if (r2 > Rsphere2)// # no => continue
				continue;
			/* Distance from the ray extreme points */
			d = alpha-dm2;
			if (d<0.)
				d=-d;
			d = d-dm2;
			if (d>0.) {
				r2 += d*d;
				/* Is this ray still able to intersect the block ? */
				if (r2 > Rsphere2)
					continue;
			}

			/* Pile initialisation with current block */
			head=0;
			igrid_pile[head] = igrid_root;
			icell_pile[head] = icell_root;
			level_pile[head] = level;
			while (head !=-1)
			{
				igrid = igrid_pile[head];
				icell = icell_pile[head];
				ilevel = level_pile[head];
				head -= 1;
#ifdef DEBUG_PRINT
                p1 = cell_centers[igrid * twotondimMndim + icell * ndim] - 0.710938;
                if (p1 < 0.0) p1 = -p1;
                p2 = cell_centers[igrid * twotondimMndim + icell * ndim + 1] - 0.164062;
                if (p2 < 0.0) p2 = -p2;
                p3 = cell_centers[igrid * twotondimMndim + icell * ndim + 2] - 0.132812;
                if (p3 < 0.0) p3 = -p3;

		        if (p1 < 1.0e-3 && p2 < 1.0e-3 && p3 < 1.0e-3)// && ilevel == 6) // && active_mask[igrid])
		        {
		            printf("Found grid (active : %d, level : %d)!\n", active_mask[igrid], level);
		        }
#endif
				/* Current cell size */
				dx = 1.0/(1<<ilevel);
				/* Cube optical depth computation */
				cube_depth = 0.0;
				intersect = 0;
				alpha = 0.0;
				alpha_min = 0.0;
				alpha_max = 0.0;
				ind_son1 = 0;
				ind_son2 = 0;
				fair_intersection = 0;
				for (idim2 = 0 ; idim2 < ndim ; idim2++) // Loop over the dimensions to search any ray-cell intersection
				{
				    jray2 = iray * ndim + idim2;
					if(ray_vects[jray2] == 0.0)
					{
						/* ray is parallel to the face along dimension idim2. There won't be an intersection. */
						continue;
					}

					/* All ray-cell intersections found => break loop */
					if (intersect == 2 && fair_intersection) break;

					/* Loop over (-,+) faces along dimension idim2 */
					for (sgn_face = -1 ; sgn_face < 2 ; sgn_face+=2) //for sgn_face in [-1,1]
					{
#ifdef DEBUG_PRINT
				    	printf("idim2 = %d ; sgn_face = %d\n", idim2, sgn_face);
#endif
						/* Coordinate parameter of the intersection of the ray with the face plane */
						alpha = ((cell_centers[igrid * twotondimMndim + icell * ndim + idim2] + sgn_face * dx / 2.0)
							 - ray_orig[jray2])/ray_vects[jray2];
						/* Is this intersection within the face ? */
						v=0.0;
						vl=0.0;
						ind_son = (twotondim-1) - (((-sgn_face+1)>>1)<<idim2);
						for (idim = 0 ; idim < ndim ; idim++){
							if (idim != idim2)
							{
							    jray = iray * ndim + idim;
								l = (ray_orig[jray] + alpha * ray_vects[jray]
								     - cell_centers[igrid * twotondimMndim + icell * ndim + idim])/dx;

								if ((l>vl) && (ray_vects[jray] == 0.0))
									vl=l;
								if (l<0.0)
								{
									l=-l;
									ind_son -= (1<<idim);
								}
								if (l>v)
									v=l;
							}
						}

						// Ensure alpha value is in [0.0, ray_front_clip] range.
						if (alpha < 0.0)
						    alpha = 0.0;
						else if (alpha > ray_front_clip)
							alpha = ray_front_clip;

#ifdef DEBUG_PRINT
				    	printf("alpha = %g, v = %g, vl = %g\n", alpha, v, vl);
#endif
						if ((v <= 0.5) && (vl < 0.5)){
							/* The intersection is in the face of the cube */
							intersect +=1;
#ifdef DEBUG_PRINT
							printf("ray_front_clip = %g\n", ray_front_clip);
							printf("intersect found : %d. Face %d (%d)\n", intersect, idim2, sgn_face);
//							for (idim = 0 ; idim < ndim ; idim++)
//							    printf("block_center[%d] = %g\n", idim, block_center[idim]);
							printf("block_size = %g\n", dx);

//							for (idim = 0 ; idim < ndim ; idim++)
//							    printf("intersection[%d] = %g\n", idim, ray_orig[iray * ndim + idim] + alpha *
//							                                            ray_vects[iray * ndim + idim]);
#endif
							if (intersect == 1) // First intersection found
							{
								alpha_min = alpha;
#ifdef DEBUG_PRINT
								printf("[intersect = 1] alpha_min = %g\n", alpha_min);
#endif
								ind_son1 = ind_son;

								continue;
							}
							else if (intersect == 2) // Second intersection found
							{
#ifdef DEBUG_PRINT
								printf("[intersect = 2] alpha = %g\n", alpha);
#endif
								alpha_max = alpha;

								if (alpha_max == alpha_min) // second intersection found is a cell corner => discard
								{
								    intersect -= 1;
								    continue;
								}
								else if (alpha_max < alpha_min) /* ind_son is the entry point son index
								                                   ind_son1 is the exit point son index */
								{
								    /* Make sure ind_son1 = entry point son index
								                 ind_son2 = exit point son index */
									ind_son2 = ind_son1;
									ind_son1 = ind_son;

									// Swap alpha_min and alpha_max to ensure alpha_max > alpha_min
									alpha_mem = alpha_min;
									alpha_min = alpha_max;
									alpha_max = alpha_mem;
								}
								else // alpha_max > alpha_min && ind_son is the exit point son index
									ind_son2 = ind_son;

								cube_depth = alpha_max - alpha_min;
#ifdef DEBUG_PRINT
								printf("[intersect = 2] cube_depth : %g\n", cube_depth);
#endif
							}
							else /* More intersections found => check we did not found the same cell corners in the
							        first two intersections */
							{
#ifdef DEBUG_PRINT
							    printf("Another intersection found : alpha = %g\n", alpha);
#endif

							    if (alpha > alpha_max) // depth = [alpha_min; alpha] && exit point son index => ind_son
							    {
							        alpha_max = alpha;
    							    cube_depth = alpha_max - alpha_min;
    							    ind_son2 = ind_son;
    							}
    							else if (alpha_min > alpha) // depth = [alpha; alpha_max] && entry point son index => ind_son
    							{
    							    alpha_min = alpha;
    							    cube_depth = alpha_max - alpha_min;
    							    ind_son1 = ind_son;
    							}
    							else
    							{
    							    intersect -= 1;
    							    continue;
    							}

#ifdef DEBUG_PRINT
    							printf("----> Another intersection found : cube_depth = %g\n", cube_depth);
#endif

    							intersect -= 1;
							}

							if (cube_depth >= dx)
							{
								fair_intersection = 1;
								break;
							}
						}
					}
				}
				/* Is the current cell actually intersected by the ray ? */
				if (cube_depth > 0.0)
				{
#ifdef DEBUG_PRINT
				    printf("Found ray-cell intersection : cube_depth : %g\n                              alpha_min : %g\n                              alpha_max : %g\n", cube_depth, alpha_min, alpha_max);
#endif
					/* The ray cuts the block */
					cell_index = igrid * twotondim + icell;
					igrid_son = sons[cell_index];

#ifdef DEBUG_PRINT
    		        printf("igrid = %d ; active_mask = %d ; igrid_son = %d\n", igrid, active_mask[igrid], igrid_son);
    		        printf("cell_center = [%g, %g, %g] ; ilevel = %d\n", cell_centers[igrid * twotondimMndim + icell * ndim],
    		                                                             cell_centers[igrid * twotondimMndim + icell * ndim + 1],
    		                                                             cell_centers[igrid * twotondimMndim + icell * ndim + 2],
    		                                                             ilevel);
#endif

					dist0 = depth_max - alpha_max;
					dist1 = depth_max - alpha_min;

					/* Test whether the cube solid angle is smaller than the ray solid angle */
					if ((igrid_son < 0) || (dist0*dist0*dS_sr[iray] >= dx*dx))
					{/* We're done searching: the subcell is the one we want */
						// We skip inactive cells
						if (active_mask[igrid] == 0)
						{
#ifdef DEBUG_PRINT
						    printf("        !!! not active cell!!!\n");
#endif
			                continue;
			            }

#ifdef DEBUG_PRINT
			            printf("                                 -> Integr. length : %g\n", dist1-dist0);
                        printf("%g ; %g ; %g  #_0_#\n", alpha_min, alpha_max, cube_depth);
#endif

						/* Add cell content to the maps */
						for (iscal = 0 ; iscal < nscal ; iscal++){
							if (use_dx)
								v = ilevel * 1.0;
							else
								v = scal_data[igrid * twotondimMnscal + icell * nscal + iscal];
							if (max_op){
								if (v > maps[iray * nscal + iscal])
									maps[iray * nscal + iscal] = v;
							}
							else
								maps[iray * nscal + iscal] += v * (dist1-dist0) * (dist1*dist1 + dist0*dist0)*dS_sr[iray]/2.0;
						}

						/* Add cube depth to the computed ray length */
						ray_l_cumul[iray] += cube_depth;
					}
					else{ /* Go down the AMR tree */
#ifdef DEBUG_PRINT
					    printf("   => ##### Go down the AMR tree ########\n");
#endif
						/* ######################################################################################### *
						 * #                      Add the NECESSARY son cells to the pile                          # *
						 * ######################################################################################### */
						if (ind_son1 != ind_son2)
						{
							for (idim = 0 ; idim < ndim ; idim++){
								ind_entry = (ind_son1>>idim) & 1;
								ind_exit = (ind_son2>>idim) & 1;
								if (ind_entry==ind_exit){
									dim_constr[idim] = 1;
									dim_side[idim] = ind_entry;
								}
								else
									dim_constr[idim] = 0;
							}

							/* Add the POSSIBLY intersecting son cells different from the entry/exit point son cells */
							for (ind_son = 0 ; ind_son < twotondim ; ind_son++){ // # All the son cells
								add = 1;
								/* Different from entry/exit point son cells */
								if ((ind_son == ind_son1)||(ind_son == ind_son2))
									add = 0;
								else
								{
									for (idim = 0 ; idim < ndim ; idim++){
										if ((dim_constr[idim])&&(dim_side[idim] != ((ind_son>>idim) & 1))){
											add = 0;
											break;
										}
									}
								}

								if (add)
								{
									head += 1;
									level_pile[head] = ilevel+1;
									igrid_pile[head] = igrid_son;
									icell_pile[head] = ind_son;
								}
							}

							/* If different from the son cell of the entry point, add the exit point son cell */
							head += 1;
							level_pile[head] = ilevel+1;
							igrid_pile[head] = igrid_son;
							icell_pile[head] = ind_son2;
						}

						/* Add the son cell of the entry point of the ray in the father grid */
						head+=1;
						level_pile[head] = ilevel+1;
						igrid_pile[head] = igrid_son;
						icell_pile[head] = ind_son1;
					}
				}
			}
		}
	}

	cleanup:

	FREE_IF_NOTNULL(igrid_pile);
	FREE_IF_NOTNULL(icell_pile);
	FREE_IF_NOTNULL(level_pile);
	FREE_IF_NOTNULL(dim_constr);
	FREE_IF_NOTNULL(dim_side);
	FREE_IF_NOTNULL(block_center);

	if (ierr == RAYTRACE_MEMORY_ERROR)
	{
		PyErr_NoMemory();
		//Py_RETURN_NONE;
		return NULL;
	}

	Py_RETURN_NONE;
}


static PyMethodDef RayTraceMethods[] = {
    { "raytrace_amr", raytrace_amr , METH_VARARGS, "Raytracing through AMR data" },
    { "raytrace_amr_spherical", raytrace_amr_spherical , METH_VARARGS, "Raytracing through AMR data (spherical projection)" },

    { NULL, NULL, 0, NULL } /* Sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,   /*m_base*/
    "_raytrace", /*m_name*/
    module_docstring,/*m_doc*/
    -1,/*m_size*/
    RayTraceMethods/*m_methods*/};


/* Python module initialisation function */
PyMODINIT_FUNC
PyInit__raytrace(void)
{
	PyObject *m = PyModule_Create(&moduledef);

	/* IMPORTANT: this must be called */
	import_array();

	if (m == NULL)
	    return NULL;

    /* Add new raytracing error type in the module */
    _RayTracingError = PyErr_NewException("_raytrace.RayTracingError", NULL, NULL);
    Py_XINCREF(_RayTracingError);
    if (PyModule_AddObject(m, "RayTracingError", _RayTracingError) < 0) {
        Py_XDECREF(_RayTracingError);
        Py_CLEAR(_RayTracingError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
