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
# Module:  mdl_kernel_regression_in_shape_space.R
# =============================================================================
#
# This module contains the methods that compute the shape regression.
# There are two main methods. Both methods require the same arguments and 
# perform the same task, the only difference being the speed.
# The two main methods are the following:
# * compute_kernel_regression_of_varset: 
#     This method uses only R code.
# * accel_compute_kernel_regression_of_varset: 
#     This method uses R code and also some C code to accelerate the most 
#     costly part of the computation. The C code is embedded by using the
#     Call interface. The C code must be compiled with the following command:
#     R CMD SHLIB. Then, the library generated (file with extension ".so") 
#     must be placed in the same directory as the R code. Then, the library
#     must be loaded with "dyn.load".
#
# This module also contain some auxiliary routines used inside these two
# methods.
#
# =============================================================================


# =============================================================================
compute_kernel_regression_of_varset = function( dataset, numSteps, varset ) {
#
# Purpose:
#   It computes a regression in the shape space of the set of variables 
#   in "varset" by using the "dataset" received. The iterative process is 
#   limited to a maximum of "numSteps" iterations.
#
# Method's arguments:
#   dataset:   Dataset to be employed in the regression.
#   numSteps:  Maximum number of steps (iterations) in the iterative process.
#   varset:    Set of variables that define the shape to be regressed.
#
  # Libraries to be used.
  library( mvtnorm )

  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # Some initializations.
  k              = dataset$k
  n              = dataset$n
  m              = dataset$m
  allCoordinates = dataset$allCoordinates
  allPreshapes   = dataset$allPreshapes
  allVariables   = dataset$allVariables
  allSizes       = dataset$allSizes

  # Check that dimensions of "varset" and "allVariables" match.
  if( length( varset ) != dim( allVariables )[ 2 ] ) {
    cat( "ERROR in 'compute_kernel_regression_of_varset': " )
    cat( "Dimensions of varset and allVariables do not match\n" )
    stop()
  }

  # Set hFactor for the initial regression process.
  hFactor = 0.50
  #### hFactor = 0.25

  # Print initial message before starting the computational process.
  cat( "Running 'compute_kernel_regression_of_varset' with arguments:\n" )
  cat( "  hFactor   = ", hFactor, "\n" )
  cat( "  k         = ", k, "\n" )
  cat( "  n         = ", n, "\n" )
  cat( "  m         = ", m, "\n" )
  cat( "  numSteps  = ", numSteps, "\n" )
  cat( "  varset    = \n" )
  print( varset )

  # --------------
  # Compute size.
  # --------------

  # Set kernel to be employed.

  #### # Spheric Gaussian Kernel.
  #### h = hFactor * apply( allVariables, 2, sd, na.rm = TRUE )
  #### dif = allVariables
  #### for( i in seq( 1, n ) ) {
  ####   dif[ i, ] = ( dif[ i, ] - varset ) / h
  #### }
  #### kernelVar = apply( dif, 1, dmvnorm )

  # Ellipsoid Gaussian Kernel.
  dif = allVariables
  for( i in seq( 1, n ) ) {
    dif[ i, ] = ( dif[ i, ] - varset )
  }
  mc = hFactor * cov( dif )
  kernelVar = apply( dif, 1, dmvnorm, sigma = mc )
  #### print( dif )
  #### print( kernelVar )

  # Compute denominator.
  denom = sum( kernelVar )

  # Check if denom is zero.
  if( denom == 0.0 ) {
    cat( "\n\n" )
    cat( "ERROR in 'compute_kernel_regression_of_varset': " )
    cat( "The prediction could not be completed\n" )
    cat( "since the denominator (variable denom) is zero.\n" )
    cat( "NULL is returned.\n\n" )
    return( NULL )
  }

  # Compute s0.
  s0 = sum( kernelVar * allSizes ) / denom

  # ---------------
  # Compute shape.
  # ---------------

  epsilon = 1.0 / n

  # Initialize q0 to any preshape (in this case, the first one).
  q0 = allPreshapes[ , , 1 ]

  #
  # Main loop.
  #
  cat( "Starting main loop...\n" )
  printConvergenceOutput = F
  tol    = 1.0e-4
  relErr = tol
  step   = 1
  while( ( step <= numSteps )&&( tol <= relErr ) )  {
    if( printConvergenceOutput ) {
      cat( "Step: ", step, "  ", proc.time()[ "elapsed" ], "s. " )
    }

    # Compute alpha.
    alpha = rep( 0, m*k-m )
    for( i in seq( 1, n ) ) {
      logba = compute_log_of_shape( allPreshapes[ , , i ], q0 ) *
              kernelVar[ i ]
      alpha = alpha + logba
    }
    alpha = ( epsilon / denom ) * alpha

    # Compute convergence.
    if( step > 1 ) {
      absErr = norm( alpha - previousAlpha, type = "F" )
      relErr = absErr / norm( alpha, type = "F" )
    } else {
      absErr = tol
      relErr = tol
    }

    # Compute the new q.
    q1 = compute_exp_of_shape( alpha, q0 )

    # Save current alpha and q.
    q0 = q1
    previousAlpha = alpha

    # Print output on convergence, if required.
    if( printConvergenceOutput ) {
      cat( "  Abs err:  ", format( absErr, scientific = TRUE ),
           "  Rel err:  ", format( relErr, scientific = TRUE ), "\n" )
    }

    step = step + 1
  }
  cat( "The iterative process required ", step - 1, " iterations.\n" )

  # ------------------------
  # Combine size and shape.
  # ------------------------

  h <- my_new_defh( k - 1 )
  q0 <- t( h ) %*% q0

  y0 = s0 * q0

  #### # Show computed object after the combination process.
  #### # This step requires the "shapes" library.
  #### plot3d( y0,
  ####         aspect = F, axes = F, xlab = "", ylab = "", zlab = "",
  ####         size = 2, type = "p" )
  #### fileName = paste( "object_after_combination", ".pdf", sep = "" )
  #### rgl.postscript( fileName, "pdf" )
  #### cat( "\n" )

  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'compute_kernel_regression_of_varset': ",
       t2 - t1, "\n" )
  cat( "\n" )

  # Return results.
  return( y0 )
}

# =============================================================================
accel_compute_kernel_regression_of_varset = function( dataset, numSteps, 
                                                      varset ) {
#
# Purpose:
#   It computes a regression in the shape space of the set of variables 
#   in "varset" by using the "dataset" received. The iterative process is 
#   limited to a maximum of "numSteps" iterations.
#   This method has been accelerated by performing the most costly parts of the
#   computation in C with BLAS libraries.
# Restriction:
#   To use this method, the C code must be compiled and loaded with "dyn.load".
#
# Method's arguments:
#   dataset:   Dataset to be employed in the regression.
#   numSteps:  Maximum number of steps (iterations) in the iterative process.
#   varset:    Set of variables that define the shape to be regressed.
#
  # Libraries to be used.
  library( mvtnorm )

  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # Some initializations.
  k              = dataset$k
  n              = dataset$n
  m              = dataset$m
  allCoordinates = dataset$allCoordinates
  allPreshapes   = dataset$allPreshapes
  allVariables   = dataset$allVariables
  allSizes       = dataset$allSizes

  # Check that dimensions of "varset" and "allVariables" match.
  if( length( varset ) != dim( allVariables )[ 2 ] ) {
    cat( "ERROR in 'accel_compute_kernel_regression_of_varset': " )
    cat( "Dimensions of varset and allVariables do not match\n" )
    stop()
  }

  # Set hFactor for the initial regression process.
  hFactor = 0.50

  # Print initial message before starting the computational process.
  cat( "Running 'accel_compute_kernel_regression_of_varset' with arguments:\n" )
  cat( "  hFactor   = ", hFactor, "\n" )
  cat( "  k         = ", k, "\n" )
  cat( "  n         = ", n, "\n" )
  cat( "  m         = ", m, "\n" )
  cat( "  numSteps  = ", numSteps, "\n" )
  cat( "  varset    = \n" )
  print( varset )

  # --------------
  # Compute size.
  # --------------

  # Set kernel to be employed.

  #### # Spheric Gaussian Kernel.
  #### h = hFactor * apply( allVariables, 2, sd, na.rm = TRUE )
  #### dif = allVariables
  #### for( i in seq( 1, n ) ) {
  ####   dif[ i, ] = ( dif[ i, ] - varset ) / h
  #### }
  #### kernelVar = apply( dif, 1, dmvnorm )

  # Ellipsoid Gaussian Kernel.
  dif = allVariables
  for( i in seq( 1, n ) ) {
    dif[ i, ] = ( dif[ i, ] - varset )
  }
  mc = hFactor * cov( dif )
  kernelVar = apply( dif, 1, dmvnorm, sigma = mc )
  #### print( dif )
  #### print( kernelVar )

  # Compute denominator.
  denom = sum( kernelVar )

  # Check if denom is zero.
  if( denom == 0.0 ) {
    cat( "\n\n" )
    cat( "ERROR in 'accel_compute_kernel_regression_of_varset': " )
    cat( "The prediction could not be completed\n" )
    cat( "since the denominator (variable denom) is zero.\n" )
    cat( "NULL is returned.\n\n" )
    return( NULL )
  }

  # Compute s0.
  s0 = sum( kernelVar * allSizes ) / denom

  # ---------------
  # Compute shape.
  # ---------------

  epsilon = 1.0 / n

  # Initialize q0 to any preshape (in this case, the first one).
  q0 = allPreshapes[ , , 1 ]

  #
  # Main loop.
  #
  cat( "Starting main loop...\n" )
  printConvergenceOutput = F
  tol = 1.0e-4
  # Perform iterative process with C code.
  q0_as_vector = as.vector( q0 )
  q0_as_vector = call_stub_of_main_iterative_process(
                     n,
                     numSteps,
                     denom,
                     printConvergenceOutput,
                     tol,
                     allPreshapes,
                     q0_as_vector,
                     kernelVar )
  # Convert plain vector into a matrix.
  q0 = matrix( q0_as_vector, nrow = k - 1 )

  # ------------------------
  # Combine size and shape.
  # ------------------------

  h <- my_new_defh( k - 1 )
  q0 <- t( h ) %*% q0

  y0 = s0 * q0

  #### # Show computed object after the combination process.
  #### # This step requires the "shapes" library.
  #### plot3d( y0,
  ####         aspect = F, axes = F, xlab = "", ylab = "", zlab = "",
  ####         size = 2, type = "p" )
  #### fileName = paste( "object_after_combination", ".pdf", sep = "" )
  #### rgl.postscript( fileName, "pdf" )
  #### cat( "\n" )

  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'accel_compute_kernel_regression_of_varset': ",
       t2 - t1, "\n" )
  cat( "\n" )

  # Return results.
  return( y0 )
}

# =============================================================================
compute_log_of_shape = function( a, b ) {
  gamma  = as.vector( b )
  alpha  = as.vector( a )
  result = alpha - gamma %*% ( t( gamma ) %*% alpha )
  return( result )
}

# =============================================================================
compute_log_of_shape_2 = function( a, b ) {
  gamma = as.vector( b )
  x     = t( b ) %*% a
  xsvd  = svd( x )
  v     = xsvd$v
  u     = xsvd$u
  tt    = v %*% t( u )
  alpha = as.vector( a %*% tt )
  result = alpha - gamma %*% ( t( gamma ) %*% alpha )
  return( result )
}

# =============================================================================
compute_exp_of_shape = function( logba, b ) {
  k = dim( b )[ 1 ] + 1
  gamma  = as.vector( b )
  inside = as.vector( sqrt( 1.0 - t( logba ) %*% logba ) ) * gamma + logba
  result = matrix( inside, nrow = k - 1 )
  return( result )
}

# =============================================================================
check_exp_log = function() {
  # Test if log and exp of shapes works.
  a = my_new_preshape( allCoordinates[ , , 1 ] )
  b = my_new_preshape( allCoordinates[ , , 2 ] )

  print( dim( a ) )
  print( dim( b ) )
  logba = compute_log_of_shape( a, b )
  z = compute_exp_of_shape( logba, b )

  cat( "Dimension of z is:  ", dim( z ), "\n" )
  cat( "max of max is:         ", max( max( z - a ) ), "\n" )
  cat( "norm of Z is:          ", norm( z, type = "F" ), "\n" )
  cat( "norm of a is:          ", norm( a, type = "F" ), "\n" )
  cat( "norm of dif is:        ", norm( z - a, type = "F" ), "\n" )
  cat( "norm_dif / norm_z is:  ", 
       norm( z - a, type = "F" ) / norm( z, type = "F" ), "\n" )

  df1 = data.frame( a = a, z = z )
  print( df1[ 1:10, ] )
}

# =============================================================================
compute_riemannian_distance = function( z, w ) {
#
# Added a tryCatch because sometimes this code aborts with the following
# error message:
# "Error in eigen(t(z) %*% w) : infinite or missing values in 'x'"
#
  dist = tryCatch( 
           {
             m <- ncol( z )
             Q <- t( z ) %*% w %*% t( w ) %*% z
             ev <- eigen( t( z ) %*% w )$values
             check <- 1
             for ( i in 1:m ) {
               check <- check * ev[ i ]
             }
             ev <- sqrt( abs( eigen( Q, symmetric = TRUE )$values ) )
             if ( Re( check ) < 0 ) {
               ev[ m ] <- -ev[ m ]
             }
             riem <- acos( min( sum( ev ), 1 ) )
             return( riem )
           },
           error = function( cond ) {
                       message( "Here's the original error message:" )
                       message( cond )
                       return( NA )
           }
         )
  return( dist )
}

# ============================================================================
call_stub_of_main_iterative_process = function( n, numSteps, denom, 
    printConvergenceOutput, tol,
    allPreshapes, q0_as_vector, kernelVar ) {
#
# It calls the C code by using the ".Call" interface.
#
  intPrintConvergenceOutput = as.integer( printConvergenceOutput )

  result = .Call( "compute_shape_in_C", 
                  n,
                  numSteps,
                  denom,
                  intPrintConvergenceOutput,
                  tol,
                  allPreshapes,
                  q0_as_vector,
                  kernelVar )
  #### print( result )
  return( result )
}

