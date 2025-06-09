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

#define OCTREE_UTILS_NO_ERROR      0
#define OCTREE_UTILS_MEMORY_ERROR  1

/* Static pointer to a octree-utils error object */
static PyObject *_OctreeUtilsError;

/* Module docstring */
static char module_docstring[] = "This module provides RAMSES octree sampling functionality.";

static void sample(double * grid_centers, int * son_indices, double * scalars, int nscalars, double * point,
                   int max_search_level, int ndim, int twotondim, int *ilevel, int *igrid, int *icell,
                   double * extracted_data, int ngrids)
{
	int idim, ibit, son_ind, iscal;
	*igrid = 0;
	// Search down levels
	for (*ilevel = 1 ; *ilevel <= max_search_level ; (*ilevel)++){
		// Find the subcell
		*icell = 0;
		for (idim = 0 ; idim < ndim ; idim++){
			ibit = point[idim] > grid_centers[(*igrid)*ndim + idim];
			*icell += (ibit << idim);
		}
		son_ind = son_indices[(*igrid)*twotondim + *icell];
		// Is the subcell a leaf cell?
		if ((son_ind < 0) || (*ilevel == max_search_level)){
			// We're done searching: the subcell is the one we want
			break;
		}
		else{
			 //Go down the tree
			*igrid = son_ind;
		}
	}
	for (iscal = 0 ; iscal < nscalars ; iscal++){
		extracted_data[iscal] = scalars[iscal*twotondim*ngrids + (*igrid)*twotondim + (*icell)];
	}
}


/* ------------------------------------------------------------------------------------------------------------------ *
 * - sample_octree_dataset() routine to sample AMR data fields a given point coordinates in a Ramses octree dataset - *
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject *
sample_octree_dataset(PyObject *self, PyObject *args) {
	int ierr = OCTREE_UTILS_NO_ERROR, read_levelmax;
	int ipoint, idim, level, igrid, icell, ilevel, ivar, iscal, i_iter, max_search_level;
	int cell_level, cell_igrid, cell_icell;
	int twotondim, twotondimMndim, ndim, nvars, npoints, nscalars, ngrids;
	int interpolation, add_level, add_cell_center;
	PyArrayObject *extracted_data_n, *points_n, *scalars_n, *grid_centers_n, *cell_centers_n, *son_indices_n;
	double *extracted_data, *points, *scalars, *grid_centers, *cell_centers;
	int *son_indices;
	double cell_size, x=0.5, y=0.5, z=0.5, *cc, xbar, ybar, zbar;
	double x_y_z=0., xbar_y_z=0., xbar_ybar_z=0., x_ybar_z=0., x_y_zbar=0., xbar_y_zbar=0., xbar_ybar_zbar=0., x_ybar_zbar=0.;
	double x_y=0., xbar_y=0., x_ybar=0., xbar_ybar=0.;

	double *x0y0z0 = NULL;
	double *x0y0z1 = NULL;
	double *x0y1z0 = NULL;
	double *x0y1z1 = NULL;
	double *x1y0z0 = NULL;
	double *x1y0z1 = NULL;
	double *x1y1z0 = NULL;
	double *x1y1z1 = NULL;
    double *point = NULL;
    double *real_point = NULL;

	if(!PyArg_ParseTuple(args, "OiOiOOOiiiiii", &points_n, &npoints, &scalars_n, &nscalars, &grid_centers_n,
	                                            &cell_centers_n, &son_indices_n, &ngrids, &ndim, &read_levelmax,
	                                            &interpolation, &add_level, &add_cell_center))
		return NULL;

	/* Fetch the numpy array data C pointers */
	points = PyArray_DATA(points_n);
	scalars = PyArray_DATA(scalars_n);
	grid_centers = PyArray_DATA(grid_centers_n);
	cell_centers = PyArray_DATA(cell_centers_n);
	son_indices = PyArray_DATA(son_indices_n);

	/* Get the integer parameters */
	twotondim = 1 << ndim;
	twotondimMndim = twotondim * ndim;
	nvars = nscalars;
	if (add_level) nvars ++;
	if (add_cell_center) nvars += ndim;

    /* Create a 2D Numpy array to store the sampled data */
	npy_intp dims[2] = {npoints, nvars};
    extracted_data_n = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    extracted_data = PyArray_DATA(extracted_data_n);

	x0y0z0 = malloc(nscalars * sizeof(double));
	x0y0z1 = malloc(nscalars * sizeof(double));
	x0y1z0 = malloc(nscalars * sizeof(double));
	x0y1z1 = malloc(nscalars * sizeof(double));
	x1y0z0 = malloc(nscalars * sizeof(double));
	x1y0z1 = malloc(nscalars * sizeof(double));
	x1y1z0 = malloc(nscalars * sizeof(double));
	x1y1z1 = malloc(nscalars * sizeof(double));

	point = malloc(ndim * sizeof(double));
	real_point = malloc(ndim * sizeof(double));

	if (x0y0z0 == NULL || x0y0z1 == NULL || x0y1z0 == NULL || x0y1z1 == NULL ||
	    x1y0z0 == NULL || x1y0z1 == NULL || x1y1z0 == NULL || x1y1z1 == NULL ||
	    point == NULL || real_point == NULL)
	{
		ierr = OCTREE_UTILS_MEMORY_ERROR;
		goto cleanup;
	}

	// Loop over points
	for (ipoint = 0 ; ipoint < npoints ; ipoint++){
		for (idim = 0 ; idim < ndim ; idim++){
			point[idim] = points[ipoint*ndim + idim];
			real_point[idim] = point[idim];
		}
		sample(grid_centers, son_indices, scalars, nscalars, point, read_levelmax, ndim, twotondim,
		       &cell_level, &cell_igrid, &cell_icell, x0y0z0, ngrids);

        level = cell_level;
        igrid = cell_igrid;
        icell = cell_icell;

		if (interpolation)
		{
			max_search_level = level;
			/* loop with 2 iterations :
			 * normally all neighbour cells has the same level and we skip the second iteration but if a lower level
			 * is found : we go into the second iteration with a reduced levelmax. There might be at most 2 levels of
			 * difference (moving in diagonal by a lenght < cell_size ) so we need at worst 3 iterations. */
			for (i_iter = 0 ; i_iter < 3 ; i_iter++)
			{
				// we compute the cell_size for this level :
				cell_size = 1;
				for (ilevel = 0 ; ilevel < level; ilevel++)
					cell_size *= .5;

				/* then we correct the difference between cell centered values and cell corner values : we just have to
				 * substract cell_size / 2 to the real (original ) point coordinates */
				for (idim = 0 ; idim < ndim ; idim++)
					point[idim] = real_point[idim] - cell_size * .5;

				sample(grid_centers, son_indices, scalars, nscalars, point, max_search_level, ndim, twotondim, &level,
				       &igrid, &icell, x0y0z0, ngrids);

                // Cell center coordinates
                cc = &cell_centers[igrid*twotondimMndim + icell*ndim];

				if (level < max_search_level)
				{
					// A lower level cell has been found, so go to the next iteration using this new levelmax
					max_search_level = level;
					continue;
				}
				// we compute first interpolation factors:
				x = (point[0] - cc[0])/cell_size + .5;
				xbar = 1.0 -x;
				y = (point[1] - cc[1])/cell_size + .5;
				ybar = 1.0 - y;
				if (ndim ==3)
				{
				    z = (point[2] - cc[2])/cell_size + .5;
				    zbar = 1.0 - z;
				    x_y_z = x * y * z;
				    xbar_y_z = xbar * y * z;
				    xbar_ybar_z = xbar * ybar * z;
				    x_ybar_z = x * ybar * z;
				    x_y_zbar = x * y * zbar;
				    xbar_y_zbar = xbar * y * zbar;
				    xbar_ybar_zbar = xbar * ybar * zbar;
				    x_ybar_zbar = x * ybar * zbar;
				}
				else
				{
				    x_y = x * y;
				    xbar_y    = xbar * y;
				    x_ybar    = x * ybar;
				    xbar_ybar = xbar * ybar;
				}


				// then we can use again the local variable "point"
				if (ndim ==3)
				{
					point[0] = cc[0];
					point[1] = cc[1];
					point[2] = cc[2] + cell_size;

					sample(grid_centers, son_indices, scalars, nscalars, point, max_search_level, ndim, twotondim,
					       &level, &igrid, &icell, x0y0z1, ngrids);

					if (level < max_search_level){
						// A lower level neighbor cell has been found, so go to the next iteration using this new levelmax
						max_search_level = level;
						continue;
					}
				}
				//point[0] = cc[0];
				point[1] = cc[1] + cell_size;
				if (ndim ==3) point[2] = cc[2];

				sample(grid_centers, son_indices, scalars, nscalars, point, max_search_level, ndim, twotondim,
				       &level, &igrid, &icell, x0y1z0, ngrids);

				if (level < max_search_level)
				{
					// A lower level neighbor cell has been found, so go to the next iteration using this new levelmax
					max_search_level = level;
					continue;
				}
				if (ndim ==3)
				{
					//point[0] = cc[0];
					//point[1] = cc[1] + cell_size;
					point[2] = cc[2] + cell_size;

					sample(grid_centers, son_indices, scalars, nscalars, point, max_search_level, ndim, twotondim,
					       &level, &igrid, &icell, x0y1z1, ngrids);

					if (level < max_search_level)
					{
						// A lower level neighbor cell has been found, so go to the next iteration using this new levelmax
						max_search_level = level;
						continue;
					}
				}
				point[0] = cc[0] + cell_size;
				point[1] = cc[1];
				if (ndim ==3) point[2] = cc[2];

				sample(grid_centers, son_indices, scalars, nscalars, point, max_search_level, ndim, twotondim,
				       &level, &igrid, &icell, x1y0z0, ngrids);

				if (level < max_search_level)
				{
					// A lower level neighbor cell has been found, so go to the next iteration using this new levelmax
					max_search_level = level;
					continue;
				}
				if (ndim ==3)
				{
					//point[0] = cc[0] + cell_size;
					//point[1] = cc[1];
					point[2] = cc[2] + cell_size;

					sample(grid_centers, son_indices, scalars, nscalars, point, max_search_level, ndim, twotondim,
					       &level, &igrid, &icell, x1y0z1, ngrids);

					if (level < max_search_level)
					{
						// A lower level neighbor cell has been found, so go to the next iteration using this new levelmax
						max_search_level = level;
						continue;
					}
				}
				//point[0] = cc[0] + cell_size;
				point[1] = cc[1] + cell_size;
				if (ndim ==3) point[2] = cc[2];

				sample(grid_centers, son_indices, scalars, nscalars, point, max_search_level, ndim, twotondim,
				       &level, &igrid, &icell, x1y1z0, ngrids);

				if (level < max_search_level)
				{
					// A lower level neighbor cell has been found, so go to the next iteration using this new levelmax
					max_search_level = level;
					continue;
				}
				if (ndim ==3)
				{
					//point[0] = cc[0] + cell_size;
					//point[1] = cc[1] + cell_size;
					point[2] = cc[2] + cell_size;

					sample(grid_centers, son_indices, scalars, nscalars, point, max_search_level, ndim, twotondim,
					       &level, &igrid, &icell, x1y1z1, ngrids);

					if (level < max_search_level)
					{
						// A lower level neighbor cell has been found, so go to the next iteration using this new levelmax
						max_search_level = level;
						continue;
					}
				}
				/* No lower level neighbor cell has been found, so we don't need to go into the next iteration of the
				 * interpolate loop */
				break;
			}

			for (iscal = 0 ; iscal < nscalars ; iscal++)
			{
				if (ndim == 3)
					extracted_data[ipoint*nvars + iscal] = xbar_ybar_zbar * x0y0z0[iscal] +
												           xbar_ybar_z * x0y0z1[iscal] + xbar_y_zbar * x0y1z0[iscal] +
													       xbar_y_z * x0y1z1[iscal] + x_ybar_zbar * x1y0z0[iscal] +
													       x_ybar_z * x1y0z1[iscal] + x_y_zbar * x1y1z0[iscal] +
													       x_y_z * x1y1z1[iscal];
				else
					extracted_data[ipoint*nvars + iscal] = xbar_ybar * x0y0z0[iscal] + xbar_y * x0y1z0[iscal] +
												           x_ybar * x1y0z0[iscal] + x_y * x1y1z0[iscal];
			}
		}
		else
		{
			for (iscal = 0 ; iscal < nscalars ; iscal++)
				extracted_data[ipoint*nvars + iscal] = x0y0z0[iscal];
		}

		ivar = nscalars;
		// Add level value in extracted data array
		if (add_level)
		{
		    extracted_data[ipoint*nvars + ivar] = (double)cell_level;
		    ivar++;
		}

		// Add cell center vector in extracted data array
		if (add_cell_center)
		{
		    cc = &cell_centers[cell_igrid*twotondimMndim + cell_icell*ndim];
		    for (idim = 0; idim < ndim ; idim++)
		    {
		        extracted_data[ipoint*nvars + ivar] = cc[idim];
		        ivar ++;
		    }
		}
	}

	if (ierr != OCTREE_UTILS_NO_ERROR)
	    goto cleanup;

	return Py_BuildValue("N", (PyObject*)extracted_data_n);

	cleanup:

	Py_XDECREF(extracted_data_n);

	FREE_IF_NOTNULL(x0y0z0);
	FREE_IF_NOTNULL(x0y0z1);
	FREE_IF_NOTNULL(x0y1z0);
	FREE_IF_NOTNULL(x0y1z1);
	FREE_IF_NOTNULL(x1y0z0);
	FREE_IF_NOTNULL(x1y0z1);
	FREE_IF_NOTNULL(x1y1z0);
	FREE_IF_NOTNULL(x1y1z1);
	FREE_IF_NOTNULL(point);
	FREE_IF_NOTNULL(real_point);

	if (ierr == OCTREE_UTILS_MEMORY_ERROR)
	{
		PyErr_NoMemory();
		return NULL;
	}


	Py_RETURN_NONE;
}

static PyMethodDef OctreeUtilsMethods[] = {
    { "sample_octree_dataset", sample_octree_dataset , METH_VARARGS, "point-sampling routine on a Ramses octree dataset" },

    { NULL, NULL, 0, NULL } /* Sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,   /*m_base*/
    "_octree_utils", /*m_name*/
    module_docstring,/*m_doc*/
    -1,/*m_size*/
    OctreeUtilsMethods/*m_methods*/};


/* Python module initialisation function */
PyMODINIT_FUNC
PyInit__octree_utils(void)
{
	PyObject *m = PyModule_Create(&moduledef);

	/* IMPORTANT: this must be called */
	import_array();

	if (m == NULL)
	    return NULL;

    /* Add new octree-utils error type in the module */
    _OctreeUtilsError = PyErr_NewException("_octree_utils.RamsesOctreeUtilsError", NULL, NULL);
    Py_XINCREF(_OctreeUtilsError);
    if (PyModule_AddObject(m, "RamsesOctreeUtilsError", _OctreeUtilsError) < 0) {
        Py_XDECREF(_OctreeUtilsError);
        Py_CLEAR(_OctreeUtilsError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}