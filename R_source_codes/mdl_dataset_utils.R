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
# Module:  mdl_dataset_utils.R
# =============================================================================
#
# This module contains several utility methods related to datasets.
#
# =============================================================================


# =============================================================================
print_dataset_dimensions = function( dataset ) {
#
# It prints the dimensions of the main components of the dataset.
#
  cat( "k:   ", dataset$k, "\n" )
  cat( "n:   ", dataset$n, "\n" )
  cat( "m:   ", dataset$m, "\n" )
  cat( "Dim of allCoordinates:  ", dim( dataset$allCoordinates ), "\n" )
  cat( "Dim of allPreshapes:    ", dim( dataset$allPreshapes ), "\n" )
  cat( "Dim of allVariables:    ", dim( dataset$allVariables ), "\n" )
  cat( "Length of allSizes:     ", length( dataset$allSizes ), "\n" )
}

# =============================================================================
create_dataset_with_selected_objects = function( dataset, vecIndices ) {
#
# It builds and returns a new dataset with a smaller number of objects in the 
# received dataset by keeping only those with indices in "vecIndices". 
# If the indices in "vecIndices" are positive, those objects are included.
# If the indices in "vecIndices" are negative, those objects are excluded.
#
  cat( "Dimensions of input dataset: \n" )
  print_dataset_dimensions( dataset )

  # Copy the received dataset into the new dataset.
  newDataset = dataset

  # Keep only the selected elements.
  newDataset$allCoordinates = newDataset$allCoordinates[ , , vecIndices ]
  newDataset$allPreshapes   = newDataset$allPreshapes  [ , , vecIndices ]
  newDataset$allVariables   = newDataset$allVariables  [ vecIndices, ]
  newDataset$allSizes       = newDataset$allSizes      [ vecIndices ]
  newDataset$k              = newDataset$k
  newDataset$n              = length( newDataset$allSizes )
  newDataset$m              = newDataset$m

  cat( "Dimensions of new dataset: \n" )
  print_dataset_dimensions( newDataset )

  # Return new dataset.
  return( newDataset )
}

# =============================================================================
reduce_num_objects_in_dataset = function( dataset, numObjects ) {
#
# It builds and returns a new dataset with a smaller number of objects in the 
# received dataset by keeping only the first "numObjects" objects in the 
# received dataset.
#
  cat( "Dimensions of input dataset: \n" )
  print_dataset_dimensions( dataset )

  # Create the new dataset.
  numObjects = min( dataset$n, numObjects )
  newIndices = seq( 1, numObjects )
  newDataset = create_dataset_with_selected_objects( dataset, newIndices )

  cat( "Dimensions of new dataset: \n" )
  print_dataset_dimensions( newDataset )

  # Return new dataset.
  return( newDataset )
}

# =============================================================================
reorder_dataset = function( dataset ) {
#
# It performs a random reordering (random permutation) of the objects in the 
# received dataset and returns it.
#
# This random ordering can be very useful when the dataset is ordered 
# (e.g.houses) and the number of elements to be processed is small, because 
# in that case the cross validation could fail since all elements in one 
# sample could be very similar.
#
  cat( "Dimensions of input dataset: \n" )
  print_dataset_dimensions( dataset )

  # Obtain number of elements in the dataset.
  numObjects = dataset$n

  # Compute permutation.
  permutation = sample( x = seq( 1, numObjects ), size = numObjects, 
                        replace = FALSE )

  # Copy the received dataset into the new dataset.
  newDataset = dataset

  # Permute the new dataset.
  newDataset$allCoordinates = dataset$allCoordinates[ , , permutation ]
  newDataset$allPreshapes   = dataset$allPreshapes  [ , , permutation ]
  newDataset$allVariables   = dataset$allVariables  [ permutation,  ]
  newDataset$allSizes       = dataset$allSizes      [ permutation ]

  cat( "Dimensions of new dataset: \n" )
  print_dataset_dimensions( newDataset )

  # Return new dataset.
  return( newDataset )
}

# =============================================================================
reduce_resolution_in_dataset = function( dataset, reductionFactor ) {
#
# It builds and returns a new dataset with the same number of objects, but 
# with a smaller number of landmarks (reduced resolution of the objects).
# If "reductionFactor" is 2, half of the points are removed, and so on.
# Only some "reductionFactor" are allowed: 1, 2, 3, and 4.
#
  # Check the reductionFactor argument.
  if( ( reductionFactor != 1 )&&
      ( reductionFactor != 2 )&&
      ( reductionFactor != 3 )&&
      ( reductionFactor != 4 ) ) {
    cat( "ERROR in 'reduce_resolution_in_dataset': Wrong reductionFactor.\n" )
    stop()
  }

  cat( "Dimensions of input dataset: \n" )
  print_dataset_dimensions( dataset )

  # Copy the received dataset into the new dataset.
  newDataset = dataset

  # Reduce the resolution of input dataset.
  newDataset$allCoordinates = newDataset$allCoordinates[ 
                                  seq( 1, newDataset$k, reductionFactor ), , ]
  newDataset$k = dim( newDataset$allCoordinates )[ 1 ]

  # Compute the new preshapes.
  k = newDataset$k
  n = newDataset$n
  m = newDataset$m
  cat( "Beginning computation of preshapes...\n" )
  preshapes = array( rep( 0, n * ( m * k - m ) ), dim = c( k-1, m, n ) )
  for( i in seq( 1, n ) ) {
    preshapes[ , , i ] = my_new_preshape( newDataset$allCoordinates[ , , i ] )
  }
  newDataset$allPreshapes = preshapes

  cat( "Dimensions of new dataset: \n" )
  print_dataset_dimensions( newDataset )

  # Return new dataset.
  return( newDataset )
}

# =============================================================================
extract_one_third = function( dataset, thirdToExtract ) {
#
# It builds and returns a list with two new datasets. The first one contains 
# one third of the original data, and the second one contains the other 
# two thirds (the rest) of the original data.
# The values of "thirdToExtract" can only be 1, 2, or 3. It represents the 
# third to be extracted from the original data.
# Some info is written in the "third.txt" file to keep track of the changes.
#
  # Remove file storing the third extracted.
  unlink( "third.txt" )

  # Check the third to be extracted.
  if( ( thirdToExtract < 1 )||( thirdToExtract > 3 ) ) {
    cat( "ERROR in 'extract_one_third': Wrong thirdToExtract\n" )
    stop()
  }

  # Check that there are at least 3 objects.
  if( dataset$n < 3 ) {
    cat( "ERROR in 'extract_one_third': Too few objects in the dataset.\n" )
    stop()
  }

  # Save the number of third to extract into a text file.
  cat( "Third to extract:  ", thirdToExtract, "\n",
       "Number of objects: ", dataset$n, "\n",
       file = "third.txt", append = T, sep = "" )

  # Compute the indices of the two sets.
  numObjsThird = dataset$n %/% 3

  third1 = seq( 1,                    numObjsThird )
  third2 = seq( numObjsThird + 1,     2 * numObjsThird )
  third3 = seq( 2 * numObjsThird + 1, dataset$n )

  if( thirdToExtract == 1 ) {
    indicesThird = third1
    indicesRest  = c( third2, third3 )
  } else if( thirdToExtract == 2 ) {
    indicesThird = third2
    indicesRest  = c( third1, third3 )
  } else if( thirdToExtract == 3 ) {
    indicesThird = third3
    indicesRest  = c( third1, third2 )
  }

  # Save some info about the partitioning into a text file.
  cat( "Number of objects in extracted third: ", numObjsThird, "\n",
       file = "third.txt", append = T, sep = "" )
  cat( "Number of objects in the rest:        ", length( indicesRest ), "\n",
       file = "third.txt", append = T, sep = "" )
  cat( "Indices of objects in extracted third: ", indicesThird, "\n",
       file = "third.txt", append = T )
  cat( "Indices of objects in the rest:        ", indicesRest, "\n",
       file = "third.txt", append = T )

  # Extract the set with the selected one third of elements.
  datasetThird = create_dataset_with_selected_objects( dataset, indicesThird )

  cat( "Dimenisions of dataset with one third of elements: \n" )
  print_dataset_dimensions( datasetThird )

  # Extract the set with the remaining two thirds of elements.
  datasetRest = create_dataset_with_selected_objects( dataset, indicesRest )

  cat( "Dimenisions of dataset with the rest of elements: \n" )
  print_dataset_dimensions( datasetRest )

  # Return results.
  return( list( datasetThird = datasetThird,
                datasetRest  = datasetRest ) )
}

# =============================================================================
compare_two_datasets = function( dataset1, dataset2 ) {
#
# It compares two datasets received as arguments.
# The result of the comparison is printed on the screen.
#
  if( ( dataset1$k != dataset2$k )||
      ( dataset1$n != dataset2$n )||
      ( dataset1$m != dataset2$m ) ) {

    cat( "Dimensions k,n,m do not match\n" )
    cat( "Dimensions of dataset1: \n" )
    print_dataset_dimensions( dataset1 )
    cat( "Dimensions of dataset2: \n" )
    print_dataset_dimensions( dataset2 )

  } else if( ! all( dim( dataset1$allCoordinates ) == 
                    dim( dataset2$allCoordinates ) ) ) {

    cat( "Dimensions of allCoordinates do not match\n" )
    cat( "Dimensions of dataset1: \n" )
    print_dataset_dimensions( dataset1 )
    cat( "Dimensions of dataset2: \n" )
    print_dataset_dimensions( dataset2 )

  } else if( ! all( dim( dataset1$allPreshapes ) == 
                    dim( dataset2$allPreshapes ) ) ) {

    cat( "Dimensions of allPreshapes do not match\n" )
    cat( "Dimensions of dataset1: \n" )
    print_dataset_dimensions( dataset1 )
    cat( "Dimensions of dataset2: \n" )
    print_dataset_dimensions( dataset2 )

  } else if( ! all( dim( dataset1$allVariables ) == 
                    dim( dataset2$allVariables ) ) ) {

    cat( "Dimensions of allVariables do not match\n" )
    cat( "Dimensions of dataset1: \n" )
    print_dataset_dimensions( dataset1 )
    cat( "Dimensions of dataset2: \n" )
    print_dataset_dimensions( dataset2 )

  } else if( length( dataset1$allSizes ) != length( dataset2$allSizes ) ) {

    cat( "Dimensions of allSizes do not match\n" )
    cat( "Dimensions of dataset1: \n" )
    print_dataset_dimensions( dataset1 )
    cat( "Dimensions of dataset2: \n" )
    print_dataset_dimensions( dataset2 )

  } else if( ! all( dataset1$allCoordinates == dataset2$allCoordinates ) ) {

    cat( "Contents of allCoordinates differ\n" )

  } else if( ! all( dataset1$allPreshapes == dataset2$allPreshapes ) ) {

    cat( "Contents of allPreshapes differ\n" )

  } else if( ! all( dataset1$allVariables == dataset2$allVariables ) ) {

    cat( "Contents of allVariables differ\n" )

  } else if( ! all( dataset1$allSizes == dataset2$allSizes ) ) {

    cat( "Contents of allSizes differ\n" )

  } else {

    cat( "Both datasets have the same contents\n" )

  }
}

