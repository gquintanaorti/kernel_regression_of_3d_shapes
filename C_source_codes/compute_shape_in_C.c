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
#include "compute_shape_in_C.h"


// ============================================================================
// Definition of local prototypes.

static void compute_exp( int n, double * pLogba, double * pBase, 
                double * pResult );

static void loop_for_computing_log( int nActual, int n,
                double * pAllPreshapes, double * pQ0_as_vector, 
                double * pKernelVar, double * pDotp, double * pResult );

static void init_double_vector( int n, double val, double * pX );

static void print_double_vector( char * vectorName, int n, double * pX );


// ============================================================================
SEXP compute_shape_in_C( 
         SEXP RnActual,
         SEXP RnumSteps,
         SEXP Rdenom,
         SEXP RprintConvergenceOutput,
         SEXP Rtol,
         SEXP allPreshapes, 
         SEXP q0_as_vector, 
         SEXP kernelVar ) {
//
// This is the main routine of this module.
// It computes the shape by using C code. It does exactly the same work as 
// the R code, but using C code and the BLAS library to accelerate the 
// processing.
//
  SEXP    result;
  double  * pAllPreshapes, * pQ0_as_vector, * pKernelVar, * pDotp, * pResult;
  double  * pQ0, * pQ1, * pPreviousAlpha, * pAlpha;
  double  denom, epsilon, factor, tol, absErr, relErr;
  int     nActual, numSteps, printConvergenceOutput, n, step;
  double  d_one = 1.0;
  int     i_one = 1;

  nActual                = INTEGER_VALUE( RnActual );
  numSteps               = INTEGER_VALUE( RnumSteps );
  epsilon                = 1.0 / ( ( double ) nActual );
  denom                  = NUMERIC_VALUE( Rdenom );
  printConvergenceOutput = INTEGER_VALUE( RprintConvergenceOutput );
  tol                    = NUMERIC_VALUE( Rtol );
  n                      = LENGTH( q0_as_vector );

  // Create SEXP object to store the result to be returned.
  result = PROTECT( allocVector( REALSXP, n ) );

  // Create auxiliary vectors.
  pAlpha         = Calloc( n, double );
  pPreviousAlpha = Calloc( n, double );
  pQ0            = Calloc( n, double );
  pQ1            = Calloc( n, double );
  pDotp          = Calloc( nActual, double );

  // Some initializations.
  pAllPreshapes = REAL( allPreshapes );
  pQ0_as_vector = REAL( q0_as_vector );
  pKernelVar    = REAL( kernelVar );
  pResult       = REAL( result );

  // Copy q0_as_vector into pQ0.
  dcopy_( & n, pQ0_as_vector, & i_one, pQ0, & i_one );

  //
  // Main loop.
  //
  relErr = tol;
  step   = 1;
  while( ( step <= numSteps )&&( tol <= relErr ) ) {

    init_double_vector( n, 0.0, pAlpha );

    loop_for_computing_log( nActual, n, 
        pAllPreshapes, pQ0, pKernelVar, pDotp, pAlpha );

    factor = epsilon / denom;
    dscal_( & n, & factor, pAlpha, & i_one );

    // Compute convergence.
    if( step > 1 ) {
      //// absErr = mydnrm2ofdiff( n, pAlpha, pPreviousAlpha );
      //// relErr = absErr / mydnrm2( n, pAlpha );
      // Overlap pPreviousAlpha with the difference, since it is not used 
      // any longer.
      daxpy_( & n, & d_one, pAlpha, & i_one, pPreviousAlpha, & i_one );
      absErr = dnrm2_( & n, pPreviousAlpha, & i_one );
      relErr = absErr / dnrm2_( & n, pAlpha, & i_one );
    } else {
      absErr = tol;
      relErr = tol;
    }

    // Compute the new q.
    // q1 = compute_exp_of_shape( alpha, q0_as_vector );
    compute_exp( n, pAlpha, pQ0, pQ1 );

    // Save current values of alpha and q1.
    dcopy_( & n, pQ1, & i_one, pQ0, & i_one );
    dcopy_( & n, pAlpha, & i_one, pPreviousAlpha, & i_one );

    // Show convergence.
    if( printConvergenceOutput ) {
      Rprintf( "  Abs err: %le  Rel err: %le \n", absErr, relErr );
    }

    step++;
  }
  //// Rprintf( "Executed %d iterationis in the main loop.\n", step - 1 );

  // Copy the result.
  dcopy_( & n, pQ0, & i_one, pResult, & i_one );

  // Remove auxiliary vectors.
  Free( pAlpha );
  Free( pPreviousAlpha );
  Free( pQ0 );
  Free( pQ1 );
  Free( pDotp );

  UNPROTECT( 1 );

  return result;
}

// ============================================================================
static void compute_exp( int n, double * pLogba, double * pBase, 
                double * pResult ) {
//
// It computes the exponential of a shape.
//
  double  dotp, factor;
  int     i;
  int     i_one  = 1;

  // # ========================================================================
  // compute_exp_of_shape = function( logba, base ) {
  //   m = dim( base )[ 2 ]
  //   k = dim( base )[ 1 ] + 1
  //   gamma = as.vector( base )
  //   inside = as.vector( sqrt( 1.0 - t( logba ) %*% logba ) ) * gamma + logba
  //   Z = matrix( inside, nrow = k-1 )
  //   return( Z )
  // }

  // inside = as.vector( sqrt( 1.0 - t( logba ) %*% logba ) ) * gamma + logba
  // dotp = 0.0 ;  
  // for( i = 0; i < n; i++ ) {
  //   dotp += pLogba[ i ] * pLogba[ i ];
  // }
  dotp = ddot_( & n, pLogba, & i_one, pLogba, & i_one );

  // In this case, no BLAS-1 is used because two operations are required, 
  // whereas the task can be performed with plain C with just one loop.
  // dcopy_( & n, pLogba, & i_one, pResult, & i_one );
  // daxpy_( & n, & factor, pBase, & i_one, pResult, & i_one );
  factor = sqrt( 1.0 - dotp );
  for( i = 0; i < n; i++ ) {
    pResult[ i ] = factor * pBase[ i ] + pLogba[ i ];
  }
}

// ============================================================================
static void loop_for_computing_log( int nActual, int n,
                double * pAllPreshapes, double * pQ0_as_vector, 
                double * pKernelVar, double * pDotp, double * pResult ) {
//
// It implements the loop including the computation of log of shapes, and
// then two operations after the loop.
// This code uses several BLAS operations to accelerate the processing for
// large vectors. Most of the acceleration comes from using BLAS-2 routines
// (the two matrix-vector operations).
//
  double  sum;
  int     i;
  double  d_one  = 1.0, d_zero = 0.0;
  int     i_one  = 1;

  // #
  // # Compute alpha.
  // #
  // alpha = rep( 0, m*k-m )
  // for( i in 1 : nActual ) {
  //   logba = compute_log_of_shape( allPreshapes[ , , i ], q0 ) *
  //           kernelVar[ i ]
  //   alpha = alpha + logba
  // }

  // First Matrix-Vector Multiply.
  dgemv_( "Transpose", & n, & nActual, 
          & d_one, pAllPreshapes, & n, pQ0_as_vector, & i_one, 
          & d_zero, pDotp, & i_one );

  // Second Matrix-Vector Multiply.
  dgemv_( "No transpose", & n, & nActual, 
          & d_one, pAllPreshapes, & n, pKernelVar, & i_one, 
          & d_one, pResult, & i_one );

  // Compute common factor.
  // sum = 0.0;
  // for( i = 0; i < nActual; i++ ) {
  //   sum += pDotp[ i ] * pKernelVar[ i ];
  // }
  sum = ddot_( & nActual, pDotp, & i_one, pKernelVar, & i_one );

  // Accumulate into result.
  // for( i = 0; i < n; i++ ) {
  //   pResult[ i ] -= sum * pQ0_as_vector[ i ];
  // }
  sum = -sum;
  daxpy_( & n, & sum, pQ0_as_vector, & i_one, pResult, & i_one );
}

// ============================================================================
static void init_double_vector( int n, double val, double * pX ) {
//
// It initializes a vector "pX" of "n" components of type double to the 
// received value "val".
//
  int  i;

  for( i = 0; i < n; i++ ) {
    pX[ i ] = val;
  }
}

// ============================================================================
static void print_double_vector( char * vectorName, int n, double * pX ) {
//
// It prints the contents of vector "pX" of "n" components of type double.
//
  int  i;

  Rprintf( "Beginning of vector %s :\n", vectorName );
  for( i = 0; i < n; i++ ) {
    Rprintf( "  pX[ %d ] = %le\n", i, pX[ i ] );
  }
  Rprintf( "End of vector.\n" );
}


