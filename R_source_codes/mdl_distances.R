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
# Module:  mdl_distances.R
# =============================================================================
#
# This module contains several methods that compute some distances between
# pairs of objects, between real objects and predicted objects, etc.
#
# =============================================================================


# =============================================================================
compute_distances_between_all_pairs_of_objects = 
    function( dataset, numObjects ) {
#
# It computes the minimum, average, and maximum distance between all pairs of
# objects in the dataset. The cost of this method is proportional to the 
# square of the number of objects since all pairs must be tested.
# Only the first numObjects objects in the dataset are used.
# The results are written into the "distances_between_pairs.txt" file.
#
# Method's arguments:
#   dataset:     Dataset containing the data.
#   numObjects:  Number of objects to be employed (-1 means all).
#

  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # Print dataset dimensions.
  print_dataset_dimensions( dataset )

  # Remove output file.
  fileName = "distances_between_all_pairs.txt"
  unlink( fileName )

  # Set variable numObjects to all objects when it is -1.
  if( numObjects == -1 ) {
    numObjects = dataset$n
  }

  # Set variable numObjects to all objects when it is larger.
  if( numObjects > dataset$n ) {
    numObjects = dataset$n
  }

  # Process every pair of objects.
  dfDist = data.frame()
  cat( "numObjects: ", numObjects, "\n" )
  for( i in seq( 1, numObjects - 1 ) ) {
    for( j in seq( i + 1, numObjects ) ) {
      # Preshapes were already computed when reading the dataset.
      #### preshape1 = my_new_preshape( dataset$allCoordinates[ , , i ] )
      #### preshape2 = my_new_preshape( dataset$allCoordinates[ , , j ] )
      preshape1 = dataset$allPreshapes[ , , i ]
      preshape2 = dataset$allPreshapes[ , , j ]
      dist = compute_riemannian_distance( preshape1, preshape2 )
      cat( "  dist between ", i, " and ", j, " is: ", dist, "\n" )

      # Add new distance to data frame.
      currRow = c( i, j, dist )
      dfDist = rbind( dfDist, currRow )
    }
  }
  colnames( dfDist ) = c( "i", "j", "dist" )
  print( dfDist )

  write.table( format( dfDist, decimal.mark = "." ),
               file = fileName,
               sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE,
               append = FALSE )

  cat( "min:  ", min(  dfDist[ , 3 ], na.rm = T ), "\n" )
  cat( "mean: ", mean( dfDist[ , 3 ], na.rm = T ), "\n" )
  cat( "max:  ", max(  dfDist[ , 3 ], na.rm = T ), "\n" )

  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'compute_distances_between_all_pairs_of_objects': ",
       t2 - t1, "\n" )
}

# =============================================================================
compare_initial_houses_and_predicted_houses_for_several_number_of_steps = 
    function( dataset ) {
#
# It computes the distances between real shapes and predicted shapes of the
# houses for several different number of steps.
# It can be used to compute the contents of table 1.
#
# Method's arguments:
#   dataset:  Dataset containing the data.
#
  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # Print dataset dimensions.
  print_dataset_dimensions( dataset )

  # Remove output file.
  fileName = "distances_between_houses_table_1.txt"
  unlink( fileName )

  # Test several number of iterations.
  for( numSteps in c( 100, 250, 500, 1000, 2000, 3000 ) ) {
    cat( "\n\n" )
    cat( "% ==============================================================\n" )
    cat( "Testing number of steps: ", numSteps, "\n" )
    cat( "% ==============================================================\n" )
    cat( "\n" )

    # Compute distances.
    dfDist = compare_initial_houses_and_predicted_houses( dataset, numSteps )

    # Write number of iterations into file.
    strLine = paste( "Testing number of steps: ", numSteps )
    write( strLine, file = fileName, append = T )

    # Write obtained data frame into file.
    strLine = paste( "i, j, k, dist" )
    write( strLine, file = fileName, append = T )
    for( i in seq( 1, nrow( dfDist ) ) ) {
      strLine = paste( dfDist[ i, 1 ], ", ", 
                       dfDist[ i, 2 ], ", ", 
                       dfDist[ i, 3 ], ", ", 
                       dfDist[ i, 4 ], sep = "" )
      write( strLine, file = fileName, append = T )
    }
    # Write summary into file.
    strLine = paste( "Min: ", min( dfDist[ , 4 ], na.rm = T ) )
    write( strLine, file = fileName, append = T )
    strLine = paste( "Med: ", mean( dfDist[ , 4 ], na.rm = T ) )
    write( strLine, file = fileName, append = T )
    strLine = paste( "Max: ", max( dfDist[ , 4 ], na.rm = T ) )
    write( strLine, file = fileName, append = T )
    strLine = paste( "" )
    write( strLine, file = fileName, append = T )

    cat( "\n\n" )
  }
  
  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'compare_distances_between_houses': ", 
       t2 - t1, "\n" )
}

# =============================================================================
compare_initial_houses_and_predicted_houses = function( dataset, numSteps ) {
#
# It generates eight perfect houses and computes the predicted houses 
# obtained in the regression process by using the dataset.
#
# Method's arguments:
#   dataset:   Dataset containing the data.
#   numSteps:  Number of steps in the iterative regression process.
#
  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # Print dataset dimensions.
  print_dataset_dimensions( dataset )

  # Some initializations.
  # IMPORTANT: This number must match the number of landmarks per edge 
  # employed when generating the dataset. Otherwise, an error might be raised 
  # at runtime.
  numLandmarksPerEdge = 3

  # Main loop.
  dfDist = data.frame()
  # for( i in c( 10 ) ) {
  #   for( j in c( 10 ) ) {
  #     for( k in c( 20 ) ) {
  for( i in c( 10, 20 ) ) {
    for( j in c( 10, 20 ) ) {
      for( k in c( 10, 20 ) ) {
        varset1 = c( i, j, k )

        cat( "\n" )
        cat( "Processing varset: ", varset1, "\n" )
  
        # Define the house for varset1 without any error.
        noErrorHouse = generate_a_house(
                           generate_model_house( numLandmarksPerEdge ),
                           numLandmarksPerEdge, varset1, 0.0 )
  
        # Predict the house for the current varset.
        # Several methods can be used here by commenting out some lines.

        ## qq0 = compute_regression_of_varset_with_tangent_space( dataset, 
        ##           varset1 )
        ## qq0 = compute_kernel_regression_of_varset( dataset, numSteps, 
        ##           varset1 )
        qq0 = accel_compute_kernel_regression_of_varset( dataset, numSteps, 
                  varset1 )
  
        # Compute both preshapes.
        predictedPreshape = my_new_preshape( qq0 )
        noErrorPreshape   = my_new_preshape( noErrorHouse )

        # Compute distance between original preshape and predicted preshape.
        dist = compute_riemannian_distance( predictedPreshape, noErrorPreshape )
        cat( "dist: ", dist, "\n" )

        #### # Show initial house.
        #### display_house( noErrorHouse )
        #### rgl.postscript( "no_error_house.pdf", "pdf" )
        #### cat( "\n" )
  
        #### # Show predicted house. 
        #### display_house( qq0 )
        #### rgl.postscript( "predicted_house.pdf", "pdf" )

        # Add new distance to data frame.
        currRow = c( i, j, k, dist )
        dfDist = rbind( dfDist, currRow )
      }
    }
  }
  colnames( dfDist ) = c( "i", "j", "k", "dist" )
  print( dfDist )

  cat( "min:  ", min(  dfDist[ , 4 ], na.rm = T ), "\n" )
  cat( "mean: ", mean( dfDist[ , 4 ], na.rm = T ), "\n" )
  cat( "max:  ", max(  dfDist[ , 4 ], na.rm = T ), "\n" )
  
  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'compare_initial_house_and_predicted_house': ",
       t2 - t1, "\n" )

  # Return data frame with results.
  return( dfDist )
}

# =============================================================================
avg_num_of_children_at_less_than_2h = function( dataset ) {
#
# It computes the average number of observations at a distance less than or 
# equal to 2h of each child.
# It can be used to generate Table 7.
#
# Method's arguments:
#   dataset:  Dataset containing the data.
#
  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # IMPORTANT: Next value must match the one used in the kernel computation 
  # routine.
  hFactor = 0.50

  # Some initializations.
  vl = vector( "numeric", dataset$n )
  sx = cov( dataset$allVariables )
  df = data.frame()

  # Process every children.
  for( i in seq( 1, dataset$n ) ) {
    d2 = mahalanobis( dataset$allVariables,
                      array( as.numeric( dataset$allVariables[ i, ] ) ),
                      sx )
    d2 = sqrt( d2 )
    numObjects = sum( d2 <= 2.0 * hFactor, na.rm = TRUE )
    vl[ i ] = numObjects
  }
  print( vl )
  avg = mean( vl )
  cat( "Avg of the number of children for h: ", hFactor, " is: ", avg, "\n" )

  currRow = c( hFactor, avg )
  df = rbind( df, currRow )

  names( df ) = c( "hFactor", "avg" )
  print( df )

  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'avg_num_of_children_at_less_than_2h': ", 
       t2 - t1, "\n" )
}

