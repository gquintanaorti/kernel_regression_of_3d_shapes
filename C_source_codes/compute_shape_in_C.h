//
// ============================================================================
// Authors
// ============================================================================
// Gregorio Quintana-Orti
//   Depto. de Ingenieria y Ciencia de Computadores,
//   Universitat Jaume I,
//   12.071 Castellon, Spain
// Amelia Simo
//   Depto. de Matematicas,
//   Universitat Jaume I,
//   12.071 Castellon, Spain
//
// ============================================================================
// Copyright
// ============================================================================
// Copyright (C) 2018,
//   Universitat Jaume I.
//
// ============================================================================
// Disclaimer
// ============================================================================
// This code is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED.
//
//
#include <R.h>
#include <Rdefines.h>
#include <R_ext/BLAS.h>
#include <math.h>

SEXP compute_shape_in_C( 
         SEXP RnActual,
         SEXP RnumSteps,
         SEXP Rdenom,
         SEXP RprintConvergenceOutput,
         SEXP Rtol,
         SEXP allPreshapes, 
         SEXP q0_as_vector, 
         SEXP kernelVar );

