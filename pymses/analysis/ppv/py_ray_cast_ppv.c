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

#define RAY_CAST_PPV_NO_ERROR      0
#define RAY_CAST_PPV_MEMORY_ERROR  1

/* Static pointer to a PPV ray-casting error object */
static PyObject *_RayCastPPVError;

static char module_docstring[] = "This module provides PPV datacube ray-casting functionality.";

/* ------------------------------------------------------------------------------------------------------------------ *
 * ---------------- ray_cast_ppv() routine to compute PPV datacubes from Ramses hydrodynamical AMR data ------------- *
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject *
ray_cast_ppv(PyObject *self, PyObject *args) {
	int ierr = RAY_CAST_PPV_NO_ERROR;
	int rlev, nblocks, nrays;
	double vmin, vmax;
	npy_intp twotondim, nx, ny, nv;
	PyArrayObject *ppv_cube_n, *density_n, *velocity_n, *ray_l_cumul_n, *ray_orig_n, *ray_vects_n, *ray_lengths_n;
	PyArrayObject *cell_centers_n, *sons_n, *igrids_n, *icells_n, *active_mask_n, *grid_levels_n;
	double *ppv_cube, *density, *velocity, *ray_l_cumul, *ray_orig, *ray_vects, *ray_depth_max, *cell_centers;
	int *sons, *igrids, *grid_levels;
	int *active_mask;
	char *icells;

	int twotondimMndim, ndim = 3;
	int igrid_root, level, iblock, idim, idim2, iray, iv, head;
	int igrid, icell, igrid_son, ilevel;
	int ind_son, intersect, sgn_face, block_already_done;
	int ind_son1, ind_son2, ind_entry, ind_exit, add, iblock2;
	char icell_root;
	double alpha, alpha_save, dx, dv, cube_depth, r2, l, d, v, rho, vl;
	double Rsphere, Rsphere2;
	double sqrt3 = sqrt(3.0);
	double dm2, depth_max;

	int* igrid_pile = NULL;
	int* icell_pile = NULL;
	int* level_pile = NULL;
	int* dim_constr = NULL;
	int* dim_side = NULL;
	double* block_center = NULL;

	if(!PyArg_ParseTuple(args, "OOOOOOOddOOiOOiOO", &ppv_cube_n, &density_n, &velocity_n, &ray_l_cumul_n, &ray_orig_n,
	                                                 &ray_vects_n, &ray_lengths_n, &vmin, &vmax, &cell_centers_n, &sons_n,
	                                                 &nblocks, &igrids_n, &icells_n, &rlev, &active_mask_n, &grid_levels_n))
		return NULL;

	dv = vmax - vmin;

	/* Get the array dimensions */
	twotondim = PyArray_DIM(density_n, 1);
	nx = PyArray_DIM(ppv_cube_n, 0);
	ny = PyArray_DIM(ppv_cube_n, 1);
	nv = PyArray_DIM(ppv_cube_n, 2);
	nrays = nx * ny;
	twotondimMndim = twotondim * ndim;

	/* Fetch the numpy array data C pointers */
	ppv_cube = PyArray_DATA(ppv_cube_n);
	density = PyArray_DATA(density_n);
	velocity = PyArray_DATA(velocity_n);
	ray_l_cumul = PyArray_DATA(ray_l_cumul_n);
	ray_orig = PyArray_DATA(ray_orig_n);
	ray_vects = PyArray_DATA(ray_vects_n);
	ray_depth_max = PyArray_DATA(ray_lengths_n);
	cell_centers = PyArray_DATA(cell_centers_n);
	sons = PyArray_DATA(sons_n);
	igrids = PyArray_DATA(igrids_n);
	icells = PyArray_DATA(icells_n);
	grid_levels = PyArray_DATA(grid_levels_n);
	active_mask = PyArray_DATA(active_mask_n);

	igrid_pile = malloc((twotondim-1) * rlev * sizeof(int));
	icell_pile = malloc((twotondim-1) * rlev * sizeof(int));
	level_pile = malloc((twotondim-1) * rlev * sizeof(int));
	dim_constr = malloc(ndim * sizeof(int));
	dim_side = malloc(ndim * sizeof(int));
	block_center = malloc(ndim * sizeof(double));

	if (igrid_pile == NULL || icell_pile == NULL || level_pile == NULL || dim_constr == NULL || dim_side == NULL ||
		block_center == NULL)
	{
		ierr = RAY_CAST_PPV_MEMORY_ERROR;
		goto cleanup;
	}

	for (iblock = 0 ; iblock < nblocks ; iblock++){
		igrid_root = igrids[iblock];
		icell_root = icells[iblock];
		block_already_done = 0;
		for (iblock2 = 0 ; iblock2 < iblock ; iblock2++){
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
			block_center[idim] = cell_centers[igrid_root * twotondimMndim + icell_root * ndim+ idim];
		}

		/* ################# *
		 * # Ray iteration # *
		 * ################# */
		for (iray = 0 ; iray < nrays ; iray++){
			depth_max = ray_depth_max[iray];
			dm2 = depth_max/2.0;
			/* Coordinate parameter of the point of the ray nearest to the block center */
			alpha = 0.0;
			for (idim = 0 ; idim < ndim ; idim++){
				alpha += (block_center[idim]-ray_orig[iray * ndim + idim])*ray_vects[iray * ndim + idim];
			}

			/* Square distances of the ray from the block center */
			r2 = 0.0;
			for (idim = 0 ; idim < ndim ; idim++){
				v = ray_orig[iray * ndim + idim] + alpha * ray_vects[iray * ndim + idim] - block_center[idim];
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
			while (head !=-1) {
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
					if(ray_vects[iray * ndim + idim2] == 0.0)
					{
						/* ray is parallel to the face along dimension idim. There might not be an intersection. */
						continue;
					}

					/* All ray-cell intersections found => break loop */
					if (intersect==2) break;
					for (sgn_face = -1 ; sgn_face < 2 ; sgn_face+=2){ /* for sgn_face in [-1,1]: */
						/* Coordinate parameter of the intersection of the ray with the face plane */
						alpha = ((cell_centers[igrid * twotondimMndim + icell * ndim + idim2] + sgn_face * dx / 2.0)
							 - ray_orig[iray * ndim + idim2])/ray_vects[iray * ndim + idim2];
						/* Is this intersection within the face ? */
						v=0.0;
						vl=0.0;
						ind_son = (twotondim-1) - (((-sgn_face+1)>>1)<<idim2);
						for (idim = 0 ; idim < ndim ; idim++){
							if (idim != idim2) {
								l = (ray_orig[iray * ndim + idim] + alpha * ray_vects[iray * ndim + idim]
								     - cell_centers[igrid * twotondimMndim + icell * ndim + idim])/dx;
								if ((l>vl) && (ray_vects[iray * ndim + idim] == 0.0))
									vl=l;
								if (l<0.0){
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
				if (cube_depth > 0.0){
					/* The ray cuts the block */
					igrid_son = sons[igrid * twotondim + icell];
					if ((igrid_son < 0) || (ilevel >= rlev))
					{/* We're done searching: the subcell is the one we want */
						// We skip inactive cells
						if (active_mask[igrid] == 0)
								continue;

						/* Add cell to PPV histogram */
						v = velocity[igrid * twotondim + icell];
						if (v >= vmin && v <= vmax)
						{
							iv = nv * (v - vmin) / dv;
							if (iv == nv) iv = nv - 1;

							rho = density[igrid * twotondim + icell];
							ray_l_cumul[iray] += cube_depth;

							ppv_cube[iray * nv + iv] += rho * cube_depth;
						}
					}
					else{ /* Go down the AMR tree */
						/* ######################################################################################### *
						 * #                      Add the NECESSARY son cells to the pile                          # *
						 * ######################################################################################### */
						if (ind_son1 != ind_son2){
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
							for (ind_son = 0 ; ind_son < twotondim ; ind_son++) // # All the son cells
							{
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

								if (add){
									head += 1;
									level_pile[head] = ilevel+1;
									igrid_pile[head] = igrid_son;
									icell_pile[head] = ind_son;
								}
							}

							/* If different from the son cell of the entry point, add the exit point son cell */
							head+=1;
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

	if (ierr == RAY_CAST_PPV_MEMORY_ERROR)
	{
		PyErr_NoMemory();
		//Py_RETURN_NONE;
		return NULL;
	}

	Py_RETURN_NONE;
}

static PyMethodDef RayCastPPVMethods[] = {
    { "ray_cast_ppv", ray_cast_ppv , METH_VARARGS, "Ray-cast a PPV datacube" },

    { NULL, NULL, 0, NULL } /* Sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,   /*m_base*/
    "_ray_cast_ppv", /*m_name*/
    module_docstring,/*m_doc*/
    -1,/*m_size*/
    RayCastPPVMethods/*m_methods*/};


/* Python module initialisation function */
PyMODINIT_FUNC
PyInit__ray_cast_ppv(void)
{
	PyObject *m = PyModule_Create(&moduledef);

	/* IMPORTANT: this must be called */
	import_array();

	if (m == NULL)
	    return NULL;

    /* Add new PPV ray-casting error type in the module */
    _RayCastPPVError = PyErr_NewException("_ray_cast_ppv.RayCastPPVError", NULL, NULL);
    Py_XINCREF(_RayCastPPVError);
    if (PyModule_AddObject(m, "RayCastPPVError", _RayCastPPVError) < 0) {
        Py_XDECREF(_RayCastPPVError);
        Py_CLEAR(_RayCastPPVError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}