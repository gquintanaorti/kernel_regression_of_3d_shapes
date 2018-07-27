#
# =============================================================================
# Authors
# =============================================================================
# Gregorio Quintana-Orti
#   Depto. de Ingenieria y Ciencia de Computadores,
#   Universitat Jaume I,
#   12.071 Castellon, Spain
# Amelia Simo
#   Depto. de Matematicas,
#   Universitat Jaume I,
#   12.071 Castellon, Spain
#
# =============================================================================
# Copyright
# =============================================================================
# Copyright (C) 2018,
#   Universitat Jaume I.
#
# =============================================================================
# Disclaimer
# =============================================================================
# This code is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED.
#


# =============================================================================
# Module:  mdl_tests.R
# =============================================================================
#
# This module contains a test of the code.
#
# =============================================================================


# =============================================================================
test_accelerated_version = function() {
#
# It tests the accelerated method by comparing against the one based on R code.
#
  # Remove all the loaded data.
  rm( list = ls_no_fun() ) 

  ds = read_dataset( "houses1" )

  varset1 = c( 12, 22, 22 )
  #### varset1 = ds$allVariables[ 2, ] + ds$allVariables[ 4, ]
  cat( "Varset to use in the regression: ", varset1, "\n\n" )

  #
  # Computing shape regression with usual R code.
  # =============================================
  #
  cat( "Computing shape regression with usual R code: \n" )

  shape1 = compute_kernel_regression_of_varset( ds, 1000, varset1 )

  #
  # Computing shape regression with accelerated code.
  # =================================================
  #
  cat( "Computing shape regression with accelerated code: \n" )

  dyn.load( "compute_shape_in_C.so" )
  shape2 = accel_compute_kernel_regression_of_varset( ds, 1000, varset1 ) 

  cat( "Comparing both results (the two shapes)...\n" )
  absErr = max( shape1 - shape2 )
  cat( "  Max( shape1 - shape2 ):                      ", absErr, "\n" )
  relErr = max( shape1 - shape2 ) / max( shape1 )
  cat( "  max( shape1 - shape2 ) / max( shape1 ):      ", relErr, "\n" )
  dist1 = compute_riemannian_distance( my_new_preshape( shape1 ), 
                                       my_new_preshape( shape2 ) )
  cat( "  Riemannian distance( preshape1, preshape2 ): ", dist1, "\n" )
  cat( "End of test_accelerated_version\n\n" )
}

