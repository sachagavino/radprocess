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

#define POINT_UTILS_NO_ERROR       0
#define POINT_UTILS_MEMORY_ERROR   1
#define POINT_NOT_A_LIST_ERROR     2
#define POINT_EMPTY_LIST_ERROR     3
#define POINT_INVALID_ARRAY_ERROR  4
#define POINT_INVALID_1DARR_ERROR  5
#define POINT_INVALID_DBLARR_ERROR 6

/* Static pointer to a point-utils error object */
static PyObject *_PointUtilsError;

/* Module docstring */
static char module_docstring[] = "This module provides meshgrid functionality.";

/* ------------------------------------------------------------------------------------------------------------------ *
 * - meshgrid() routine to build an array of point coordinates fromed by the cartesian product of the arrays listed - *
 * - the input array.                                                                                               - *
 * -                                                                                                                - *
 * - The points are listed in C (row-major) order.                                                                  - *
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject *
meshgrid(PyObject *self, PyObject *args) {
	int ierr = POINT_UTILS_NO_ERROR;
	int ndim, idim, icoord;
	PyObject *coords_by_axis, *list_item;
	PyArrayObject *arr_item, *points_n = NULL;
	double *points, *darr;
	long npoints = 1, ipoint, iperiod, irepeat;
	long *shape = NULL;
	long *nperiods = NULL;
	long *nrepeats = NULL;

	if(!PyArg_ParseTuple(args, "O", &coords_by_axis))
		return NULL;

    // Check coords_by_axis is a valid list
    if (! PyList_Check(coords_by_axis))
    {
        ierr = POINT_NOT_A_LIST_ERROR;
        goto cleanup;
    }

	/* Extract the list of coordinate arrays */
	ndim = (int)PyList_Size(coords_by_axis);
	if(ndim <= 0)
	{
	    ierr = POINT_EMPTY_LIST_ERROR;
	    goto cleanup;
	}

	/* Allocate shape/stride arrays */
	shape = malloc(ndim * sizeof(long));
    nperiods = malloc(ndim * sizeof(long));
    nrepeats = malloc(ndim * sizeof(long));
	if (shape == NULL ||nperiods == NULL || nrepeats == NULL)
	{
		ierr = POINT_UTILS_MEMORY_ERROR;
		goto cleanup;
	}

	/* Fill shape/stride arrays */
	for (idim = 0; idim < ndim; idim++)
	{
		list_item = PyList_GetItem(coords_by_axis, idim);
		if(!PyArray_Check(list_item))
		{
		    ierr = POINT_INVALID_ARRAY_ERROR;
		    goto cleanup;
		}
		arr_item = (PyArrayObject*)list_item;
		if (PyArray_NDIM(arr_item) != 1)
		{
		    ierr = POINT_INVALID_1DARR_ERROR;
		    goto cleanup;
		}
		if (PyArray_TYPE(arr_item) != NPY_DOUBLE)
		{
		    ierr = POINT_INVALID_DBLARR_ERROR;
		    goto cleanup;
		}
		nperiods[idim] = npoints;
		shape[idim] = (long)PyArray_DIM(arr_item, 0);
		npoints *= shape[idim];
	}
	for (idim = 0; idim < ndim; idim++)
		nrepeats[idim] = npoints / (nperiods[idim] * shape[idim]);

    /* Create a 2D Numpy array to store the point coordinates */
	npy_intp dims[2] = {npoints, ndim};
    points_n = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    points = PyArray_DATA(points_n);

	// For every coordinate...
	for (idim = 0; idim < ndim; idim++)
	{
		// Get all the values which the coord idim takes... slow, but in outer loop
		arr_item = (PyArrayObject*)PyList_GetItem(coords_by_axis, idim);
		darr = (double*)PyArray_DATA(arr_item);

		ipoint = 0;
		icoord = 0;
		for (iperiod = 0; iperiod < nperiods[idim]; iperiod++)
		{
			for (icoord = 0; icoord < shape[idim]; icoord++)
			{
				for (irepeat = 0 ; irepeat < nrepeats[idim]; irepeat++)
				{
					points[ipoint*ndim + idim] = darr[icoord];
					ipoint ++;
				}
			}
		}
    }

	if (ierr != POINT_UTILS_NO_ERROR)
	    goto cleanup;

	return Py_BuildValue("N", (PyObject*)points_n);

	cleanup:

	Py_XDECREF(points_n);

	FREE_IF_NOTNULL(shape);
	FREE_IF_NOTNULL(nperiods);
	FREE_IF_NOTNULL(nrepeats);

	if (ierr == POINT_UTILS_MEMORY_ERROR)
	{
		PyErr_NoMemory();
		return NULL;
	}
	else if(ierr == POINT_NOT_A_LIST_ERROR)
	{
		PyErr_SetString(_PointUtilsError, "'coords_by_axis' attribute is not a valid list object.");
	    return NULL;
	}
	else if(ierr == POINT_EMPTY_LIST_ERROR)
	{
		PyErr_SetString(_PointUtilsError, "'coords_by_axis' list attribute must contain at least one coordinate vector.");
	    return NULL;
	}
	else if(ierr == POINT_INVALID_ARRAY_ERROR)
	{
		PyErr_SetString(_PointUtilsError, "'coords_by_axis' list must only contain Numpy array elements.");
	    return NULL;
	}
	else if(ierr == POINT_INVALID_1DARR_ERROR)
	{
		PyErr_SetString(_PointUtilsError, "'coords_by_axis' list must only contain coordinate vector (1D Numpy array) elements.");
	    return NULL;
	}
	else if(ierr == POINT_INVALID_DBLARR_ERROR)
	{
		PyErr_SetString(_PointUtilsError, "'coords_by_axis' list must only contain (double value) coordinate vectors.");
	    return NULL;
	}

	Py_RETURN_NONE;
}

static PyMethodDef PointUtilsMethods[] = {
    { "meshgrid", meshgrid , METH_VARARGS, "n-dimensional point coordinate list computation from coordinate vector list" },

    { NULL, NULL, 0, NULL } /* Sentinel */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,   /*m_base*/
    "_point_utils", /*m_name*/
    module_docstring,/*m_doc*/
    -1,/*m_size*/
    PointUtilsMethods/*m_methods*/};


/* Python module initialisation function */
PyMODINIT_FUNC
PyInit__point_utils(void)
{
	PyObject *m = PyModule_Create(&moduledef);

	/* IMPORTANT: this must be called */
	import_array();

	if (m == NULL)
	    return NULL;

    /* Add new point-utils error type in the module */
    _PointUtilsError = PyErr_NewException("_point_utils.PointUtilsError", NULL, NULL);
    Py_XINCREF(_PointUtilsError);
    if (PyModule_AddObject(m, "PointUtilsError", _PointUtilsError) < 0) {
        Py_XDECREF(_PointUtilsError);
        Py_CLEAR(_PointUtilsError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}