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
# Module:  mdl_cross_validation.R
# =============================================================================
#
# This module performs cross validations.
#
# =============================================================================


# =============================================================================
cross_validation = function( dataset ) {
#
# It performs a cross validation of the received dataset.
# The dataset is partitioned extracting both one third and the rest to perform 
# a 3-fold cross validation.
#
# Method's arguments:
#   dataset:  Dataset to be processed.
#

  #
  # Perform a random ordering of the received dataset.
  # This random ordering can be very useful when the dataset is ordered 
  # (e.g.houses) and the number of elements to be processed is small, because 
  # in that case the cross validation could fail since all elements in one 
  # third could be similar.
  #
  # Set seed for reproducibility.
  set.seed( 11 )
  dataset = reorder_dataset( dataset )

  # Reduce the data set.
  # To perform some shorter tests, the total number of objects can be reduced.
  dataset = reduce_num_objects_in_dataset( dataset, 30 )

  # Prepare output file name.
  fileName = paste( "output_of_cross_validation_with_", 
                    dataset$n, "_objects.txt", sep = "" )

  # Remove output file.
  unlink( fileName )

  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # Process several number of iterations.
  for( numSteps in c( 100, 250, 500, 1000, 2000 ) ) {
    cat( "\n\n\n" )
    cat( "-----------------------------------------------------------------\n" )
    cat( "Processing numSteps: ", numSteps, "\n" )
    cat( "-----------------------------------------------------------------\n" )

    # Write an initial message to output file.  
    cat( "\n", file = fileName, append = T )
    cat( "------------------------------------\n", file = fileName, append = T )
    cat( "Processing numSteps: ", numSteps, "\n", file = fileName, append = T )
    cat( "------------------------------------\n", file = fileName, append = T )
    cat( "\n", file = fileName, append = T )

    # Process the three thirds.
    vecDist = rep( NA, length.out = dataset$n )
    i = 1
    for( third in seq( 1, 3 ) ) {
      cat( "\n" )
      cat( "---------------------------------------------------------------\n" )
      cat( "Processing third: ", third, "\n" )
      cat( "---------------------------------------------------------------\n" )
  
      # Particiona dataset en dos datasets: current third and the rest.
      tlist              = extract_one_third( dataset, third )
      dsetOfCurrentThird = tlist$datasetThird
      dsetOfRest         = tlist$datasetRest
   
      for( j in seq( 1, dsetOfCurrentThird$n ) ) {
        #
        # Process every object inside the current third.
        #
        # Extract variables of object "j".
        cat( "Variables of object j: ", j, 
             " in dsetOfCurrentThird: ", third, "\n" )
        objectVars = dsetOfCurrentThird$allVariables[ j, ]
        print( objectVars )
        varset = as.vector( as.matrix( objectVars ) )
  
        # Predict shape of object "j" out of its varset.
        qq = compute_kernel_regression_of_varset( dsetOfRest, numSteps, 
                 varset )
        ### qq = accel_compute_kernel_regression_of_varset( dsetOfRest, 
        ###          numSteps, varset )
  
        # Check if prediction was obtained 
        if( is.null( qq ) ) {
          # No prediction was obtained
          dist = NA
        } else { 
          # Prediction was obtained
          # Compute both preshapes.
          originalPreshape = dsetOfCurrentThird$allPreshapes[ , , j ]
          newPreshape      = my_new_preshape( qq )
    
          # Compute distance between original preshape and predicted preshape.
          dist = compute_riemannian_distance( originalPreshape, newPreshape )
        }
        cat( "dist: ", dist, "\n\n" )
    
        # Add new distance to vector.
        vecDist[ i ] = dist 
  
        i = i + 1
      }
      # Print data so far.
      cat( "Partial results including third:", third, "\n" )
      cat( "Vector of distances:\n" )
      print( vecDist )
      vmini = min(  vecDist, na.rm = T )
      vmean = mean( vecDist, na.rm = T )
      vmaxi = max(  vecDist, na.rm = T )
      cat( "Third ", third, " so-far min:  ", vmini, "\n" )
      cat( "Third ", third, " so-far mean: ", vmean, "\n" )
      cat( "Third ", third, " so-far max:  ", vmaxi, "\n" )
  
      # Save partial results into file.
      cat( "Partial results including third:", third, "\n", 
           file = fileName, append=T )
      cat( "Vector of distances:\n", file = fileName, append = T )
      cat( vecDist, "\n", file = fileName, append = T )
      cat( "Partial min:  ", vmini, "\n", file = fileName, append = T )
      cat( "Partial mean: ", vmean, "\n", file = fileName, append = T )
      cat( "Partial max:  ", vmaxi, "\n", file = fileName, append = T )
      cat( "\n", file = fileName, append = T )
    }
  
    # Print total data.
    cat( "\n" )
    cat( "Total results:\n" )
    cat( "vecDist:\n" )
    print( vecDist )
    vmini = min(  vecDist, na.rm = T )
    vmean = mean( vecDist, na.rm = T )
    vmaxi = max(  vecDist, na.rm = T )
    cat( "Total min:  ", vmini, "\n" )
    cat( "Total mean: ", vmean, "\n" )
    cat( "Total max:  ", vmaxi, "\n" )
   
    # Save final results into file.
    cat( "\n", file = fileName, append = T )
    cat( "Total results for numSteps : ", numSteps, "\n", 
         file = fileName, append = T )
    cat( "vecDist: \n", file = fileName, append = T )
    cat( vecDist, "\n", file = fileName, append = T )
    cat( "Total min:  ", vmini, "\n", file = fileName, append = T )
    cat( "Total mean: ", vmean, "\n", file = fileName, append = T )
    cat( "Total max:  ", vmaxi, "\n", file = fileName, append = T )
    cat( "\n", file = fileName, append = T )
 
  } 

  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'validacion_cruzada': ",
       t2 - t1, "\n" )
}

