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
#include "read_amr.h"
#include "read_cells.h"
#include "read_parts.h"
#include "utils.h"
#include <strings.h>

/* Static pointer to a Ramses I/O error object */
static PyObject *_ReadRamsesError;

/* Module docstring */
static char module_docstring[] = "This module provides low-level RAMSES output I/O functionalities.";
/* ----------------------------------------------------------------------------------------------------------------- */
/* ----------------------------- Shortcut routines for easy dictionary insertion ----------------------------------- */
/* ----------------------------------------------------------------------------------------------------------------- */
void dict_int(PyObject * dict, const char * key, int val) {
	PyObject *valobj;
	valobj = PyLong_FromLong(val);
	PyDict_SetItemString(dict, key, valobj);
	Py_DECREF(valobj);
}

void dict_float(PyObject * dict, const char * key, double val) {
	PyObject *valobj;
	valobj = PyFloat_FromDouble(val);
	PyDict_SetItemString(dict, key, valobj);
	Py_DECREF(valobj);
}

void dict_ndarr(PyObject * dict, const char * key, void * data, int nd, npy_intp * dims, int type) {

	/* Compute total element count */
	int size = 1;
	int i;
	for(i=0; i<nd; i++) size *= dims[i];

	/* Create empty array */
	PyObject *arrobj;
	arrobj = PyArray_SimpleNew(nd, dims, type);

	/* Copy into array buffer */
	void *npdata = PyArray_DATA((PyArrayObject*)arrobj);
	memcpy(npdata, data, PyArray_ITEMSIZE((PyArrayObject*)arrobj)*size);

	/* Add to dict */
	PyDict_SetItemString(dict, key, arrobj);
	Py_DECREF(arrobj);
}

void dict_1darr(PyObject * dict, const char * key, void * data, npy_intp n, int type) {
	dict_ndarr(dict, key, data, 1, &n, type);
}


void dict_2darr(PyObject * dict, const char * key, void * data, npy_intp nx, npy_intp ny, int type) {
	npy_intp dims[2] = {nx, ny};
	dict_ndarr(dict, key, data, 2, dims, type);
}
/* ----------------------------------------------------------------------------------------------------------------- */




/* ------------------------------------------------------------------------------------------------------------------ *
 * ------------------ read_amr() routine to import data from Ramses 'amr_XXXXX.outYYYY' F90 binary files ------------ *
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject * read_amr(PyObject *self, PyObject *args) {
	const char *filename;
	char err_string[256];
	int readlmax;
	int swap = 0;
	int ierr = READ_AMR_NO_ERROR;
	RAMSES_AmrHeader_t * hdr;
	RAMSES_AmrStruct_t * amr;
	PyObject *amr_dict, *hdr_dict;

	if(!PyArg_ParseTuple(args, "si|i", &filename, &readlmax, &swap))
		return NULL;


    /* PyErr_SetString(_ReadRamsesError, "Invalid attributes for read_amr() function"); */

	/* Read the AMR structure */
	hdr = RAMSES_AmrHeader_New();
	amr = RAMSES_AmrStruct_New();
	if (hdr == NULL || amr == NULL) {
	    ierr = READ_AMR_MEMORY_ERROR;
	    goto cleanup;
	}

	ierr = RAMSES_ReadAmr(hdr, amr, filename, readlmax, swap);

	if (ierr != READ_AMR_NO_ERROR) goto cleanup;

	/* Assemble the header dictionary */
	hdr_dict = PyDict_New();

	dict_int   (hdr_dict, "ndim"         , hdr->ndim       );
	dict_int   (hdr_dict, "ncpu"         , hdr->ncpu       );

	dict_1darr (hdr_dict, "coarse_shape" , hdr->coarse_shape, 3, NPY_INT);

	dict_int   (hdr_dict, "ncoarse"      , hdr->ncoarse    );
	dict_int   (hdr_dict, "nboundary"    , hdr->nboundary  );
	dict_int   (hdr_dict, "levelmax"     , hdr->levelmax   );
	dict_int   (hdr_dict, "ngridmax"     , hdr->ngridmax   );
	dict_int   (hdr_dict, "noutputs"     , hdr->noutputs   );

	dict_1darr (hdr_dict, "out_times"    , hdr->out_times, hdr->noutputs, NPY_DOUBLE);
	dict_1darr (hdr_dict, "out_aexps"    , hdr->out_aexps, hdr->noutputs, NPY_DOUBLE);

	dict_float (hdr_dict, "boxlen"       , hdr->boxlen     );
	dict_int   (hdr_dict, "iout"         , hdr->iout       );
	dict_float (hdr_dict, "aexp"         , hdr->aexp       );
	dict_float (hdr_dict, "time"         , hdr->time       );
	dict_int   (hdr_dict, "ngrids_lmax"  , hdr->ngrids_lmax);

	/* Assemble the AMR structure dictionary */
	amr_dict = PyDict_New();

	dict_int   (amr_dict, "ndim"              , amr->ndim       );
	dict_int   (amr_dict, "twotondim"         , amr->twotondim  );
	dict_int   (amr_dict, "ncpu"              , amr->ncpu       );
	dict_int   (amr_dict, "nboundary"         , amr->nboundary  );
	dict_int   (amr_dict, "levelmax"          , amr->levelmax   );
	dict_int   (amr_dict, "ngrids"            , amr->ngrids     );
	dict_int   (amr_dict, "ngrids_lmax"       , amr->ngrids_lmax);
	dict_int   (amr_dict, "readlmax"          , amr->readlmax   );

	dict_2darr (amr_dict, "ngridlevel"        , amr->ngridlevel        , amr->ncpu     , amr->levelmax , NPY_INT);
	dict_2darr (amr_dict, "ngridbound"        , amr->ngridbound        , amr->nboundary, amr->levelmax , NPY_INT);
	dict_1darr (amr_dict, "grid_indices"      , amr->grid_indices      , amr->ngrids                   , NPY_INT);
	dict_1darr (amr_dict, "coarse_son_indices", amr->coarse_son_indices, amr->ncoarse                  , NPY_INT);
	dict_2darr (amr_dict, "grid_centers"      , amr->grid_centers      , amr->ngrids   , amr->ndim     , NPY_DOUBLE);
	dict_2darr (amr_dict, "son_indices"       , amr->son_indices       , amr->ngrids   , amr->twotondim, NPY_INT);

	RAMSES_AmrStruct_Free(amr);
	RAMSES_AmrHeader_Free(hdr);

	return Py_BuildValue("(N,N)", hdr_dict, amr_dict);

	/* Cleanup */
	cleanup:

	RAMSES_AmrStruct_Free(amr);
	RAMSES_AmrHeader_Free(hdr);

	if (ierr == READ_AMR_MEMORY_ERROR) {
		PyErr_NoMemory();
		//Py_RETURN_NONE;
		return NULL;
	}
	else if (ierr == RAED_AMR_NGRID_MISMATCH_ERROR) {
	    sprintf(err_string, "Number of read AMR grids mismatch while reading AMR file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;

	}
	else if (ierr == FIO_OPEN_FAIL) {
		sprintf(err_string, "Failed to open F90 binary file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_CLOSE_FAIL) {
		sprintf(err_string, "Failed to close F90 binary file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_READ_RECORD_FAIL) {
		sprintf(err_string, "Failed to read F90 record from file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_SKIP_RECORD_FAIL) {
		sprintf(err_string, "Failed to skip F90 record from file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	Py_RETURN_NONE;
}

/* ------------------------------------------------------------------------------------------------------------------ *
 * -------- read_hydro_header() routine to import metadata from Ramses AMR hydro files 'hydro_XXXXX.outYYYY', ------- *
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject * read_hydro_header(PyObject *self, PyObject *args) {
	const char *filename;
	char err_string[256];
	int swap = 0;
	int ierr = READ_CELLS_NO_ERROR;
	RAMSES_CellHydroHeader_t *hdr = NULL;

	if(!PyArg_ParseTuple(args, "s|i", &filename, &swap))
		return NULL;

	/* Allocate the header C structure */
	hdr = RAMSES_CellHydro_Header_New();
	if (hdr == NULL) {
		ierr = READ_CELLS_MEMORY_ERROR;
		goto cleanup;
	}

	/* Read the Hydro file header structure */
	ierr = RAMSES_CellHydro_Header_Read(hdr, filename, swap);

	/* Assemble the header dictionary */
	PyObject *hdr_dict = PyDict_New();

	dict_int  (hdr_dict, "ndim"     , hdr->ndim     );
	dict_int  (hdr_dict, "ncpu"     , hdr->ncpu     );
	dict_int  (hdr_dict, "nvar_file", hdr->nvar_file);
	dict_int  (hdr_dict, "levelmax" , hdr->levelmax );
	dict_float(hdr_dict, "gamma"    , hdr->gamma    );

	if (ierr != READ_CELLS_NO_ERROR)
		goto cleanup;

	RAMSES_CellHydro_Header_Free(hdr);

	return Py_BuildValue("N", hdr_dict);

	/* Cleanup */
	cleanup:

	RAMSES_CellHydro_Header_Free(hdr);

	if (ierr == READ_CELLS_MEMORY_ERROR)
	{
		PyErr_NoMemory();
		//Py_RETURN_NONE;
		return NULL;
	}
	else if (ierr == FIO_OPEN_FAIL)
	{
		sprintf(err_string, "Failed to open F90 binary file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_CLOSE_FAIL)
	{
		sprintf(err_string, "Failed to close F90 binary file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_READ_RECORD_FAIL)
	{
		sprintf(err_string, "Failed to read F90 record from file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_SKIP_RECORD_FAIL)
	{
		sprintf(err_string, "Failed to skip F90 record from file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}

	Py_RETURN_NONE;
}


/* ------------------------------------------------------------------------------------------------------------------ *
 * -------- read_cells() routine to import data from Ramses AMR physical field files like 'hydro_XXXXX.outYYYY', ---- *
 * -------------------- 'grav_XXXXX.outYYYYY', 'mhd_XXXXX.outYYYY', etc. F90 binary files --------------------------- *
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject *
read_cells(PyObject *self, PyObject *args) {
	const char *filename;
	char err_string[256];
	int ngrids;
	int is_hydro_file;
	int ndim;
	int readlmax;
	int swap = 0;
	int grav_compat = 0;
	int ierr = READ_CELLS_NO_ERROR;
	RAMSES_CellData_t * cell = NULL;
	PyArrayObject *data_array = NULL;

	PyObject * ivarlist;
	if(!PyArg_ParseTuple(args, "siiiiO|ii", &filename, &is_hydro_file, &ndim, &ngrids, &readlmax, &ivarlist, &swap,
	                     &grav_compat))
		return NULL;

	/* Extract the list of ivars to read */
	Py_ssize_t nvarout = PyList_Size(ivarlist);
	if(nvarout <= 0)
	{
	    ierr = READ_CELLS_NODATA_ERROR;
	    goto cleanup;
	}

	/* Allocate cell data structure */
	cell = RAMSES_CellData_New((int)nvarout, ndim, ngrids);
	if (cell == NULL)
	{
		ierr = READ_CELLS_MEMORY_ERROR;
		goto cleanup;
	}

	/* Fill variable index array */
	Py_ssize_t ielem;
	for (ielem=0; ielem<nvarout; ielem++) {
		PyObject * item = PyList_GetItem(ivarlist, ielem);
		if(!PyLong_Check(item))
		{
		    ierr = READ_CELLS_INVALID_IDATA_ERROR;
		    goto cleanup;
		}
		long ivar = PyLong_AsLong(item);
		cell->ivar_arr[ielem] = (int)ivar;
	}

    /* Create a 3D Numpy array to store the cell data */
	npy_intp dims[3] = {nvarout, (1<<cell->ndim), ngrids};
    data_array = (PyArrayObject*)PyArray_SimpleNew(3, dims, NPY_DOUBLE);
    cell->data = (double *)PyArray_DATA(data_array);

	/* Read the cell data */
	ierr = RAMSES_CellData_Read(cell, filename, is_hydro_file, readlmax, swap, grav_compat);


	if (ierr != READ_CELLS_NO_ERROR)
	    goto cleanup;

	RAMSES_CellData_Free(cell);

	return Py_BuildValue("N", (PyObject*)data_array);

    /* Cleanup */
    cleanup:

	Py_XDECREF(data_array);

	RAMSES_CellData_Free(cell);

	if (ierr == READ_CELLS_MEMORY_ERROR)
	{
	    PyErr_NoMemory();
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if(ierr == READ_CELLS_NODATA_ERROR)
	{
		PyErr_SetString(_ReadRamsesError, "No variable index was selected.");
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == READ_CELLS_NDIM_INVALID)
	{
		PyErr_SetString(_ReadRamsesError, "'ndim' parameter is not consistent.");
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if(ierr == READ_CELLS_INVALID_IDATA_ERROR)
	{
		PyErr_SetString(_ReadRamsesError, "Variable index is not an integer value.");
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == READ_CELLS_READ_NGRID_MISS)
	{
		sprintf(err_string, "Missing grids were not read. File '%s' may be corrupted.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == READ_CELLS_READ_NGRID_OVERFLOW)
	{
		sprintf(err_string, "Trying to read to many grid values. File '%s' may be corrupted.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr >= READ_CELLS_UNAVAIL_IDATA_ERROR)
	{
		int ivar = (ierr - READ_CELLS_UNAVAIL_IDATA_ERROR);
		sprintf(err_string, "Variable of index '%d' is not available.", ivar);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_OPEN_FAIL)
	{
		sprintf(err_string, "Failed to open F90 binary file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_CLOSE_FAIL)
	{
		sprintf(err_string, "Failed to close F90 binary file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_READ_RECORD_FAIL)
	{
		sprintf(err_string, "Failed to read F90 record from file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_SKIP_RECORD_FAIL)
	{
		sprintf(err_string, "Failed to skip F90 record from file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}

	Py_RETURN_NONE;
}


/* ------------------------------------------------------------------------------------------------------------------ *
 * --------------- read_parts() routine to import data from Ramses 'part_XXXXX.outYYYY' F90 binary files ------------ *
 * ------------------------------------------------------------------------------------------------------------------ */
static PyObject *
read_parts(PyObject *self, PyObject *args) {
	const char *filename;
	char err_string[256];
	FIO_File * fiofile = NULL;
	RAMSES_PartData_t * part;
	double boxlen;
	int ierr, swap = 0, def_longint = 0;

	if(!PyArg_ParseTuple(args, "sd|ii", &filename, &boxlen, &swap, &def_longint))
		return NULL;


	/* Open file and read particle data into memory */
	fiofile = FIO_Open(filename, "r", 4, swap? FIO_SWAP_BYTES : FIO_DEFAULT);
	if (fiofile == NULL)
	{
		ierr = FIO_OPEN_FAIL;
		goto cleanup;
	}

	part = RAMSES_PartData_New();
	if (part == NULL)
	{
		ierr = READ_PARTS_MEMORY_ERROR;
		goto cleanup;
	}

	/* Read particle data from file */
	ierr = RAMSES_PartData_Read(part, boxlen, fiofile, def_longint);

	if (ierr != READ_PARTS_NO_ERROR)
	    goto cleanup;

	ierr = FIO_Close(fiofile);
	if (ierr != 0)
		goto cleanup;

	/* Header */
	PyObject *hdr_dict = PyDict_New();
	dict_int(hdr_dict, "ncpu" , part->ncpu);
	dict_int(hdr_dict, "ndim" , part->ndim);
	dict_int(hdr_dict, "npart", part->npart);
	dict_int(hdr_dict, "nstar", part->nstar);
	dict_int(hdr_dict, "nsink", part->nsink);

	/* Particle data */
	PyObject *part_dict = PyDict_New();

	if(part->pos != NULL)
		dict_2darr(part_dict, "pos"  , part->pos  , part->npart, part->ndim, NPY_DOUBLE);

	if(part->vel != NULL)
		dict_2darr(part_dict, "vel"  , part->vel  , part->npart, part->ndim, NPY_DOUBLE);

	if(part->mass != NULL)
		dict_1darr(part_dict, "mass" , part->mass , part->npart,             NPY_DOUBLE);

	if((!def_longint && part->id_32 != NULL) || (def_longint && part->id_64 != NULL))
	{
	    if (def_longint)
	        dict_1darr(part_dict, "id"   , part->id_64   , part->npart,             NPY_INT64);
        else
        	dict_1darr(part_dict, "id"   , part->id_32   , part->npart,             NPY_INT);
	}

	if(part->level != NULL)
		dict_1darr(part_dict, "level", part->level, part->npart,             NPY_INT);

	if(part->epoch != NULL)
		dict_1darr(part_dict, "epoch", part->epoch,	part->npart,             NPY_DOUBLE);

	if(part->metal != NULL)
		dict_1darr(part_dict, "metal", part->metal, part->npart,             NPY_DOUBLE);

	/* Cleanup */
	RAMSES_PartData_Free(part);

	return Py_BuildValue("(N,N)", hdr_dict, part_dict);

	cleanup:

	if (fiofile != NULL)
		FIO_Close(fiofile);

	if (ierr == READ_PARTS_MEMORY_ERROR)
	{
	    PyErr_NoMemory();
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_OPEN_FAIL)
	{
		sprintf(err_string, "Failed to open F90 binary file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_CLOSE_FAIL)
	{
		sprintf(err_string, "Failed to close F90 binary file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_READ_RECORD_FAIL)
	{
		sprintf(err_string, "Failed to read F90 record from file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}
	else if (ierr == FIO_SKIP_RECORD_FAIL)
	{
		sprintf(err_string, "Failed to skip F90 record from file '%s'.", filename);
		PyErr_SetString(_ReadRamsesError, err_string);
	    //Py_RETURN_NONE;
	    return NULL;
	}

	Py_RETURN_NONE;
}


static PyMethodDef ReadRamsesMethods[] = {
    { "read_amr",          read_amr ,          METH_VARARGS, "Read an AMR structure" },
	{ "read_hydro_header", read_hydro_header , METH_VARARGS, "Read hydro file header"},
    { "read_cells",        read_cells ,        METH_VARARGS, "Read data from AMR cells"},
    { "read_parts",        read_parts ,        METH_VARARGS, "Read RAMSES particle data"},

    { NULL, NULL, 0, NULL } /* Sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,   /*m_base*/
    "_read_ramses", /*m_name*/
    module_docstring,/*m_doc*/
    -1,/*m_size*/
    ReadRamsesMethods/*m_methods*/};


/* Python module initialisation function */
PyMODINIT_FUNC
PyInit__read_ramses(void)
{
	PyObject *m = PyModule_Create(&moduledef);

	/* IMPORTANT: this must be called */
	import_array();

	if (m == NULL)
	    return NULL;

    /* Add new Ramses I/O error type in the module */
    _ReadRamsesError = PyErr_NewException("_read_ramses.RamsesIOError", NULL, NULL);
    Py_XINCREF(_ReadRamsesError);
    if (PyModule_AddObject(m, "RamsesIOError", _ReadRamsesError) < 0) {
        Py_XDECREF(_ReadRamsesError);
        Py_CLEAR(_ReadRamsesError);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
