/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Tridiagonal matrix algorithm
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "Python.h"
#include <numpy/arrayobject.h>
#include <stdlib.h>

double *pyvector_to_Carrayptrs (PyArrayObject * arrayin);

static PyObject *
TDMASolver (PyObject * self, PyObject * args)
{
  double mc;
  double *ac, *bc, *cc, *dc, *xc;
  double *acnew, *bcnew, *ccnew, *dcnew;
  int i, npoints;

  /* Variables for calling from Python */
  PyArrayObject *a, *b, *c, *d;
  PyObject *out_array;

  /*  parse single numpy array argument with a series of double variables after it */
  if (!PyArg_ParseTuple
      (args, "O!O!O!O!", &PyArray_Type, &a, &PyArray_Type, &b,
       &PyArray_Type, &c, &PyArray_Type, &d))
    return NULL;


  /* output array to match the input arrays */
  out_array = PyArray_NewLikeArray (b, NPY_ANYORDER, NULL, 0);
  if (out_array == NULL)
    return NULL;
  acnew = calloc (sizeof (double), a->dimensions[0]);
  bcnew = calloc (sizeof (double), b->dimensions[0]);
  ccnew = calloc (sizeof (double), c->dimensions[0]);
  dcnew = calloc (sizeof (double), d->dimensions[0]);

  npoints = d->dimensions[0];


  /* initial and output arrays */
  ac = pyvector_to_Carrayptrs (a);
  bc = pyvector_to_Carrayptrs (b);
  cc = pyvector_to_Carrayptrs (c);
  dc = pyvector_to_Carrayptrs (d);
  xc = pyvector_to_Carrayptrs (out_array);

  for (i = 0; i < npoints; i++)
  {
    bcnew[i] = bc[i];
    dcnew[i] = dc[i];
  }

  for (i = 1; i < npoints; i++)
  {
    // printf ("%.1f %.1f %.1f %.1f\n", ac[i], bc[i],cc[i], dc[i]);
    mc = ac[i-1] / bcnew[i-1];
    bcnew[i] = bcnew[i] - mc * cc[i - 1]; 
    dcnew[i] = dcnew[i] - mc * dcnew[i - 1];
  }

  /* set boundary condition */
  xc[npoints - 1] = dcnew[npoints - 1] / bcnew[npoints - 1];

  /* second sweep */
  for (i = npoints - 2; i > -1; i--)
  {
    xc[i] = (dcnew[i] - cc[i] * xc[i+1]) / bcnew[i];
  }

  /*  clean up and return the result */
  free (acnew);
  free (bcnew);
  free (ccnew);
  free (dcnew);
  Py_INCREF (out_array);
  return Py_BuildValue("O", out_array);
}

/* Create 1D Carray from PyArray
   Assumes PyArray is contiguous in memory. Credit Lou Pecora     */
double *
pyvector_to_Carrayptrs (PyArrayObject * arrayin)
{
  int n;

  n = arrayin->dimensions[0];
  return (double *) arrayin->data;  /* pointer to arrayin data as double */
}


/* Python interfacing stuff */

/*  define functions in module */
static PyMethodDef TDMAMethods[] = {
  {"TDMASolver", TDMASolver, METH_VARARGS,
   "Tridiagonal matrix algorithm"},
  {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
/* module initialization */
/* Python version 3*/
static struct PyModuleDef cModTDMA = {
  PyModuleDef_HEAD_INIT,
  "TDMA_module", "Tridiagonal matrix algorithm",
  -1,
  TDMAMethods
};

PyMODINIT_FUNC
PyInit_tdma (void)
{
  PyObject *module;
  module = PyModule_Create (&cModTDMA);
  if (module == NULL)
    return NULL;
  /* IMPORTANT: this must be called */
  import_array ();
  if (PyErr_Occurred ())
    return NULL;
  return module;
}

#else
/* module initialization */
/* Python version 2 */
PyMODINIT_FUNC
inittdma (void)
{
  PyObject *module;
  module = Py_InitModule ("tdma", TDMAMethods);
  if (module == NULL)
    return;
  /* IMPORTANT: this must be called */
  import_array ();
  return;
}

#endif
