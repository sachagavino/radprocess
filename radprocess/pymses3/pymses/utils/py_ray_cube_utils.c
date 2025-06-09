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

/* Static pointer to a ray/cube utils error object */
static PyObject *_RayCubeUtilsError;

/* Module docstring */
static char module_docstring[] = "This module provides unit-sized cube coordinates clipping functionality.";

static int _ray_cube_intersect(double *ray_orig, double *ray_vect, double ray_length_max, double *cube_center,
                               double dx_cube, int ndim, double *intersect_length, double *alpha_min, double *alpha_max,
                               int *ind_son_entry, int *ind_son_exit)
{
    /**
     * Computes the number of intersections of a finite ray starting from an origin point and with a given length with a
     * 2D/3D cube, given its size and center coordinates.
     *
     * Parameters :
     * ------------
     * ray_orig: ray origin coordinate vector (size ndim)
      * ray_vect: ray direction coordinate vector (size ndim)
     * ray_length_max: maximum ray length (starting from origin point)
     * cube_center: cube center coordinate vector (size ndim)
     * dx_cube: size of the cube
     * ndim: number of dimensions (2 or 3)
     *
     * Updates :
     * ---------
     * intersect_length: length of the intersection of the finite ray with the cube (>= 0.0)
     * alpha_min:
     * alpha_max:
     * ind_son_entry:
     * ind_son_exit:
     *
     * Returns :
     * ---------
     *   * intersect = -1 if the ray does not intersect the cube :
     *           +---------+     /
     *           |         |    /
     *           |         |   /
     *           |         |  /
     *           +---------+ /
     *   * intersect = 0 if the ray is included in the cube:
     *           +---------+
     *           |    /    |
     *           |   /     |
     *           |  /      |
     *           +---------+
     *   * intersect = 1 if the ray intersects the cube through a corner or halfway through the cube:
     *      \                              \
     *       \   +---------+           +----\----+
     *        \  |         |           |     \   |
     *         \ |         |      OR   |      \  |
     *          \|         |           |         |
     *           +---------+           +---------+
     *            \
     *   * intersect = 2 if the ray pierces the cube through and through:
     *                  /
     *           +-----/---+
     *           |    /    |
     *           |   /     |
     *           |  /      |
     *           +-/-------+
     *            /
     */
    double alpha, alpha_save, dm_2, Rsphere, Rsphere2;
    int idim, idim2, sgn_face, nintersect, int_ret, twotondimMone, ind_son, ind_son_save;
	double l, d, r2, v, vl;
	double sqrt3 = sqrt(3.0);
	double dx_2 = dx_cube / 2.0;
    Rsphere = sqrt3 * dx_cube;
    Rsphere2 = Rsphere * Rsphere;
    twotondimMone = (1 << ndim) - 1;


    /* Coordinate parameter of the point of the ray nearest to the unit size cube */
	alpha = 0.0;
	for (idim = 0 ; idim < ndim ; idim++)
		alpha += (cube_center[idim] - ray_orig[idim]) * ray_vect[idim];

	/* Square distances of the ray from the unit size cube */
	r2 = 0.0;
	for (idim = 0 ; idim < ndim ; idim++)
	{
		v = ray_orig[idim] + alpha * ray_vect[idim] - cube_center[idim];
		r2 += v*v;
	}

	/* Is this ray able to intersect the cube ? */
	if (r2 > Rsphere2)// no => return -1
    {
        *intersect_length = 0.0;
        *alpha_min = 0.0;
        *alpha_max = 0.0;
        return -1;
    }

	/* Distance from the ray extreme points */
	dm_2 = ray_length_max / 2.0;
	d = alpha - dm_2;
	if (d < 0.0)
		d=-d;
	d -= dm_2;
	if (d > 0.0)
	{
		r2 += d*d;

	    /* Is this ray still able to intersect the cube ? */
    	if (r2 > Rsphere2)
    	{
    	    // no => return -1
            *intersect_length = 0.0;
            *alpha_min = 0.0;
            *alpha_max = 0.0;
            return -1;
        }
	}

	/* Cube intersection computation */
	nintersect = 0;
	alpha = 0.0;
	alpha_save = 0.0;
	ind_son = 0;
	ind_son_save = 0;

	for (idim = 0 ; idim < ndim ; idim++)
	{
	    /* All ray-cell intersections found => break loop */
		if (nintersect == 2) break;

		if(ray_vect[idim] == 0.0)
		{
			/* ray is parallel to the face along dimension idim. There won't be an intersection. */
			continue;
		}

        /* Iteration over the +/- cube face along dimension idim */
		for (sgn_face = -1 ; sgn_face < 2 ; sgn_face += 2){ /* for sgn_face in [-1,1]: */
			/* Coordinate of the intersection of the ray with the face plane */
			alpha = ((cube_center[idim] + sgn_face * dx_2) - ray_orig[idim]) / ray_vect[idim];
			/* Is this intersection within the face ? */
			v=0.0;
			vl=0.0;
			ind_son = twotondimMone - (((-sgn_face + 1) >> 1) << idim);
			for (idim2 = 0 ; idim2 < ndim ; idim2++)
			{
			    /* Intersection normalised coordinates (v, vl) absolute values in the cube face */
				if (idim != idim2)
				{
					l = (ray_orig[idim2] + alpha * ray_vect[idim2] - cube_center[idim2]) / dx_cube;
					if ((l>vl) && (ray_vect[idim2] == 0.0)) vl=l;
					if (l<0.0)
					{
					    l=-l;
					    ind_son -= (1 << idim2);
					}
		   			if (l>v) v=l;
				}
			}

            /* Is the intersection within the cube face */
			if (v <= 0.5 && vl < 0.5)
			{
			    nintersect ++;
			    if (nintersect == 1) // First intersection found
			    {
				    alpha_save = alpha;
				    ind_son_save = ind_son;
				}
    			else // Second intersection found : done
					break;
			}
		}
	}

    if (nintersect == 0) // No intersection found => return -1
    {
        *intersect_length = 0.0;
        *alpha_min = 0.0;
        *alpha_max = 0.0;
    	*ind_son_entry = 0;
	    *ind_son_exit = 0;
        return -1;
    }
	else if (nintersect == 1)
	{
        *intersect_length = 0.0;

	    /* Clip alpha_save to [0.0; ray_length_max] */
        if (alpha_save < 0.0 || alpha_save > ray_length_max)
        {
            *alpha_min = 0.0;
            *alpha_max = 0.0;
        	*ind_son_entry = 0;
    	    *ind_son_exit = 0;
            return -1;
        }

	    // Position of the single intersection with the cube along the ray (cube corner)
	    *alpha_min = alpha_save;
	    *alpha_max = alpha_save;
    	*ind_son_entry = ind_son_save;
	    *ind_son_exit = ind_son_save;
	    return 1;
	}
	else if (nintersect == 2)
	{
	    int_ret = 0;

        // Face intersection is within the finite ray ? => increase intersection counter
        if (alpha >= 0.0 && alpha <= ray_length_max) int_ret++;
        else // Face intersection is out of the finite ray => clip coordinate along the finite ray length
        {
            /* Clip alpha to [0.0; ray_length_max] */
            if (alpha < 0.0) alpha = 0.0;
            else if (alpha > ray_length_max) alpha = ray_length_max;
        }

        // Face intersection is within the finite ray ? => increase intersection counter
        if (alpha_save >= 0.0 && alpha_save <= ray_length_max) int_ret++;
        else // Face intersection is out of the finite ray => clip coordinate along the finite ray length
        {
            /* Clip alpha_save to [0.0; ray_length_max] */
            if (alpha_save < 0.0) alpha_save = 0.0;
            else if (alpha_save > ray_length_max) alpha_save = ray_length_max;
        }

	    // Sort coordinates
	    if (alpha > alpha_save)
	    {
	        *alpha_min = alpha_save;
	        *ind_son_entry = ind_son_save;
	        *alpha_max = alpha;
	        *ind_son_exit = ind_son;
	    }
	    else
	    {
	        *alpha_min = alpha;
	        *ind_son_entry = ind_son;
	        *alpha_max = alpha_save;
	        *ind_son_exit = ind_son_save;
	    }

	    // Length of the intersection with the cube
	    *intersect_length = *alpha_max - *alpha_min;

	    // Return (effective) number of intersections
	    return int_ret;
	}

	return -1;
}


static PyObject *
clip_unit_size_cube(PyObject *self, PyObject *args) {
	int ierr = RAYTRACE_NO_ERROR, needs_clipping;
	PyArrayObject *ray_orig_n, *ray_vects_n, *ray_lengths_n;
	double *ray_orig, *ray_vects, *ray_lengths;
	int ndim, nrays;
	int idim, iray, jray, ind_son_entry, ind_son_exit;
	int intersect;
	double alpha_min, alpha_max, dx=1.0;
	double *cube_center = NULL;
	double cube_depth;

	if(!PyArg_ParseTuple(args, "OOO", &ray_orig_n, &ray_vects_n, &ray_lengths_n))
		return NULL;

	/* Fetch the numpy array data C pointers */
	ray_orig = PyArray_DATA(ray_orig_n);
	ray_vects = PyArray_DATA(ray_vects_n);
	ray_lengths = PyArray_DATA(ray_lengths_n);

	/* Get the integer parameters */
	nrays = (int)PyArray_DIM(ray_orig_n, 0);
	ndim = (int)PyArray_DIM(ray_orig_n, 1);

	cube_center = malloc(ndim * sizeof(double));

	if (cube_center == NULL)
	{
		ierr = RAYTRACE_MEMORY_ERROR;
		goto cleanup;
	}

	/* ############################################################################################################## *
	 * ############################################### Ray iteration ################################################ *
	 * ############################################################################################################## */
	for (iray = 0 ; iray < nrays ; iray++)
	{
	    // Check whether the ray_origin lies outside the unit size cube
	    needs_clipping = 0;
		for (idim = 0 ; idim < ndim ; idim++)
		{
		    jray = iray * ndim + idim;
			if (ray_orig[jray] < 0.0 || ray_orig[jray] > dx)
			{
			    needs_clipping = 1;
			    break;
			}
		}
		if (! needs_clipping)
		    continue;

        // Set unit size cube center coordinates
		for (idim = 0 ; idim < ndim ; idim++)
			cube_center[idim] = 0.5;

        intersect = _ray_cube_intersect(&ray_orig[iray*ndim], &ray_vects[iray*ndim], ray_lengths[iray], cube_center, dx,
                                        ndim, &cube_depth, &alpha_min, &alpha_max, &ind_son_entry, &ind_son_exit);

        // No intersection found : set the ray_length to 0.0 and continue
        if (intersect < 0)
		{
		    ray_lengths[iray] = 0.0;
		    continue;
		}

        /* Ray is included in the unit size cube ? Should never happen since the ray_origin lays outside the unit size
         * cube => continue. */
		if (intersect == 0)
		    continue;

        // At least one intersection found...

		/* Shift the ray origin to the closest intersection of the ray with the unit size cube */
		if (alpha_min > 0.0)
		{
		    for (idim = 0 ; idim < ndim ; idim++)
		    {
    		    jray = iray * ndim + idim;
    			ray_orig[jray] += alpha_min * ray_vects[jray];
    		}
        }

        /* Update the ray_length to the computed cube_depth value */
        ray_lengths[iray] = cube_depth;
	}
	/* ############################################################################################################## */

	cleanup:

	FREE_IF_NOTNULL(cube_center);

	if (ierr == RAYTRACE_MEMORY_ERROR)
	{
		PyErr_NoMemory();
		return NULL;
	}

	Py_RETURN_NONE;
}

static PyMethodDef RayTraceCubeClipMethods[] = {
    { "clip_unit_size_cube", clip_unit_size_cube, METH_VARARGS, "Clipping ray origins to the unit size cube"},

    { NULL, NULL, 0, NULL } /* Sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,   /*m_base*/
    "_ray_cube_utils", /*m_name*/
    module_docstring,/*m_doc*/
    -1,/*m_size*/
    RayTraceCubeClipMethods/*m_methods*/};


/* Python module initialisation function */
PyMODINIT_FUNC
PyInit__ray_cube_utils(void)
{
	PyObject *m = PyModule_Create(&moduledef);

	/* IMPORTANT: this must be called */
	import_array();

	if (m == NULL)
	    return NULL;

    /* Add new ray/cube utils error type in the module */
    _RayCubeUtilsError = PyErr_NewException("_ray_cube_utils.RayCubeUtilsError", NULL, NULL);
    Py_XINCREF(_RayCubeUtilsError);
    if (PyModule_AddObject(m, "RayCubeUtilsError", _RayCubeUtilsError) < 0) {
        Py_XDECREF(_RayCubeUtilsError);
        Py_CLEAR(_RayCubeUtilsError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}