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

/* Static pointer to a Ramses I/O error object */
static PyObject *_LICError;

static char module_docstring[] = "This module provides Line Integral Convolution 2D map computation functionality.";

// Error flag definitions
#define LIC_NO_ERROR                0
#define LIC_MEMORY_ERROR            1
#define LIC_NOT_FLOAT_DTYPE_ERROR   2
#define LIC_WRONG_DIMS_ERROR        3
#define LIC_INVALID_SHAPES_ERROR    4
#define INVALID_KERNEL_ERROR        5

/* Move to the next pixel in the vector field direction (follow streamline). This function perform x, y, fx, and fy
 * in-place updates.
 *
 * Params
 * ------
 *
 *
 *   vx : float
 *     Vector x component.
 *   vy :float
 *     Vector y component.
 *   x : int
 *     Pixel x index. Updated in place.
 *   y : int
 *     Pixel y index. Updated in place.
 *   fx : float
 *     Position along x in the pixel unit square. Updated in place.
 *   fy : float
 *     Position along y in the pixel unit square. Updated in place.
 *   w : int
 *     Number of pixels along x.
 *   h : int
 *     Number of pixels along y.
 */
static void _streamline_integrate(float vx, float vy, int* x, int* y, float*fx, float*fy, int w, int h)
{
    float tx, ty;  /* Think of tx (ty) as the time it takes to reach the next pixel along the x (y) axis. */
    int zeros = 0;

    /*printf("vx = %f, vy = %f\n", vx, vy);
    printf("fx = %f, fy = %f\n", *fx, *fy);
    printf("x = %d, y = %d\n", *x, *y);
    printf("w=%d, h=%d\n", w, h);*/


    if (vx > 0)
        tx = (1- *fx) / vx;
    else if (vx < 0)
        tx = - *fx / vx;
    else
    {
        zeros ++;
        tx = 1e100;
    }

    if (vy > 0)
        ty = (1 - *fy) / vy;
    else if (vy < 0)
        ty = - *fy / vy;
    else
    {
        zeros ++;
        ty = 1e100;
    }

    if (zeros == 2)
        return;

    if (tx<ty) // We reached the next pixel along x first.
    {
        if (vx>=0)
        {
            (*x)++;
            *fx = 0;
        }
        else
        {
            (*x)--;
            *fx = 1;
        }
        *fy += tx*vy;
    }
    else // We reached the next pixel along y first.
    {
        if (vy>=0)
         {
            (*y)++;
            *fy = 0;
        }
        else
        {
            (*y)--;
            *fy = 1;
        }
        *fx += ty*vx;
    }

    /* Streamline reached image border, clip to image */
    if (*x >= w)
        *x = w-1;
    if (*x < 0)
        *x = 0;
    if (*y >= h)
        *y = h-1;
    if (*y < 0)
        *y = 0;
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * Return an image of the texture array blurred along the local vector field orientation.
 *
 * Parameters
 * ----------
 *  vel_u_n : array (ny, nx)
 *      Vector field x-axis component.
 *  vel_v_n : array (ny, nx)
 *      Vector field y-axis component.
 *  nois_texture_n : array (ny,nx)
 *    The texture image that will be distorted by the vector field. Usually, a white noise image is recommended to
 *    display the fine structure of the vector field.
 *  kernel : 1D array
 *    The convolution kernel: an array weighting the texture along the stream line. For static images, a box kernel
 *    (equal to one) of length max(nx,ny)/10 is appropriate. The kernel should be symmetric.
 * norm : int
 *    normalize the output array as to have values in the [0.0, 1.0] range. Default 1 (True).
 *
 *  Returns
 *  -------
 *  out : array(ny,nx)
 *    An image of the texture convoluted along the vector field streamlines.
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject *
texture_streamline_convol(PyObject *self, PyObject *args) {
	PyArrayObject *vel_u_n, *vel_v_n, *noise_texture_n, *kernel_n;
	float *vel_u, *vel_v, *noise_texture, *kernel, *res_map;
	npy_intp *map_size;
	PyArrayObject *res_map_n = NULL;
	char err_string[256];
	int i, j, k, x, y, nx, ny, ker_length, kmax, norm=1;
	int ierr = LIC_NO_ERROR;
    float fx, fy, tx, ty, min_pix=1e100, max_pix=-1e100;

	if(!PyArg_ParseTuple(args, "OOOO|i",&vel_u_n, &vel_v_n, &noise_texture_n, &kernel_n, &norm))
		return NULL;

    /* check input numpy array dtypes are 'float' */
    if (PyArray_TYPE(vel_u_n) != NPY_FLOAT32 || PyArray_TYPE(vel_v_n) != NPY_FLOAT32 ||
        PyArray_TYPE(noise_texture_n) != NPY_FLOAT32 || PyArray_TYPE(kernel_n) != NPY_FLOAT32)
    {
        ierr = LIC_NOT_FLOAT_DTYPE_ERROR;
        goto cleanup;
    }

    /* check input numpy arrays are 2D arrays */
    if (PyArray_NDIM(vel_u_n) != 2 || PyArray_NDIM(vel_v_n) != 2 || PyArray_NDIM(noise_texture_n) != 2)
    {
        ierr = LIC_WRONG_DIMS_ERROR;
        goto cleanup;
    }

    // Fetch map size + kernel length
    map_size = PyArray_SHAPE(noise_texture_n);
    nx = (int)map_size[0];
    ny = (int)map_size[1];
    ker_length = (int)PyArray_DIM(kernel_n, 0);

    /* Check input nnumpy arrays have identical shapes */
    if (PyArray_DIM(vel_u_n, 0) != nx || PyArray_DIM(vel_u_n, 1) != ny ||
        PyArray_DIM(vel_v_n, 0) != nx || PyArray_DIM(vel_v_n, 1) != ny)
    {
        ierr = LIC_INVALID_SHAPES_ERROR;
        goto cleanup;
    }

    /* Check that kernel function is a 1D numpy.ndarray of odd size */
    if (PyArray_NDIM(kernel_n) > 1 || ker_length % 2 != 1)
    {
        ierr == INVALID_KERNEL_ERROR;
        goto cleanup;
    }
    kmax = (ker_length - 1)/2;

	/* Fetch the numpy array data C pointers */
	vel_u = (float *)PyArray_DATA(vel_u_n);
	vel_v = (float *)PyArray_DATA(vel_v_n);
	noise_texture = (float *)PyArray_DATA(noise_texture_n);
	kernel = (float *)PyArray_DATA(kernel_n);

    // Init result map array
    res_map_n = (PyArrayObject*)PyArray_SimpleNew(2, map_size, NPY_FLOAT32);
    res_map = (float *)PyArray_DATA(res_map_n);

    for (i=0; i<nx; i++)
    {
        for (j=0; j<ny; j++)
        {
            x = i;
            y = j;
            fx = 0.5;// pixel center coordinates
            fy = 0.5;// ------------------------

            k = kmax;
            res_map[i*ny + j] = kernel[k] * noise_texture[x*ny + y];

            while (k < ker_length-1)
            {
                _streamline_integrate(vel_u[x*ny + y], vel_v[x*ny + y], &x, &y, &fx, &fy, nx, ny);
                k++;
                res_map[i*ny + j] += kernel[k] * noise_texture[x*ny + y];
            }

            x = i;
            y = j;
            fx = 0.5;// pixel center coordinates
            fy = 0.5;// ------------------------

            k = kmax;

            while (k>0)
            {
                _streamline_integrate(-vel_u[x*ny + y], -vel_v[x*ny + y], &x, &y, &fx, &fy, nx, ny);
                k--;
                res_map[i*ny + j] += kernel[k] * noise_texture[x*ny + y];
            }

            if (res_map[i*ny + j] > max_pix)
                max_pix = res_map[i*ny + j];
            if (res_map[i*ny + j] < min_pix)
                min_pix = res_map[i*ny + j];
        }
    }

    /* Final LIC image normalisation to obtain [0.0, 1.0] value range */
    if(norm) {
        for (i=0; i<nx; i++)
        {
            for (j=0; j<ny; j++)
            {
                res_map[i*ny + j] = (res_map[i*ny + j] - min_pix) / (max_pix - min_pix);
            }
        }
    }

    return Py_BuildValue("N", (PyObject*)res_map_n);

    /* Cleanup */
    cleanup:

	Py_XDECREF(res_map);

	if(ierr == LIC_NOT_FLOAT_DTYPE_ERROR)
	{
		PyErr_SetString(_LICError, "'vel_u_n', 'vel_v_n', 'noise_texture_n' and 'kernel_n' must be numpy.ndarray instances of dtype 'float'.");
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if(ierr == LIC_WRONG_DIMS_ERROR)
	{
		PyErr_SetString(_LICError, "'vel_u_n', 'vel_v_n' and 'noise_texture_n' must be 2D numpy.ndarray instances.");
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == LIC_INVALID_SHAPES_ERROR)
	{
    	PyErr_SetString(_LICError, "'vel_u_n', 'vel_v_n' and 'noise_texture_n' must have identical shapes.");
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == INVALID_KERNEL_ERROR)
	{
    	PyErr_SetString(_LICError, "'kernel_n' must be a odd length 1D numpy.ndarray instance.");
	    //Py_RETURN_NONE;
	    return NULL;
	}

	sprintf(err_string, "Error during Line integral convolution processing !");
	PyErr_SetString(_LICError, err_string);
	//Py_RETURN_NONE;
	return NULL;
}


static PyMethodDef LICMethods[] = {
    { "texture_streamline_convol", texture_streamline_convol , METH_VARARGS, "Compute streamline integral convolution 2D map" },

    { NULL, NULL, 0, NULL } /* Sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,   /*m_base*/
    "_line_integral_convolution", /*m_name*/
    module_docstring,/*m_doc*/
    -1,/*m_size*/
    LICMethods/*m_methods*/};


/* Python module initialisation function */
PyMODINIT_FUNC
PyInit__line_integral_convolution(void)
{
	PyObject *m = PyModule_Create(&moduledef);

	/* IMPORTANT: this must be called */
	import_array();

	if (m == NULL)
	    return NULL;

    /* Add new LIC error type in the module */
    _LICError = PyErr_NewException("_line_integral_convolution.LineIntegralConvolutionError", NULL, NULL);
    Py_XINCREF(_LICError);
    if (PyModule_AddObject(m, "LineIntegralConvolutionError", _LICError) < 0) {
        Py_XDECREF(_LICError);
        Py_CLEAR(_LICError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}