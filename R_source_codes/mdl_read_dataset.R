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
# Module:  mdl_read_dataset.R
# =============================================================================
#
# This module contains a function to read the dataset from two files: the 
# first one with the coordinates and the second one with the set of variables 
# (varset).
#
# =============================================================================


# =============================================================================
read_dataset = function( strTestName ) {
#
# Read the dataset contained in the two following files:
#     Data/strTestName_coor.csv
#     Data/strTestName_vars.csv
# The first file contains the coordinates of the landmarks for every object.
# The second file contains the variables for every object.
# After reading the dataset, preshapes are computed to accelerate further
# computations.
#
  # Check input argument.
  if( strTestName == "" ) {
    cat( "\n" )
    cat( "ERROR in read_dataset: test name is empty.\n" )
    cat( "Usage:  read_dataset( strTestName ).\n" )
    cat( "There must exist the following two files: \n" )
    cat( "    Data/strTestName_coor.csv\n" )
    cat( "    Data/strTestName_vars.csv\n" )
    cat( "\n" )
    stop()
  }

  # Set timer.
  t1 = proc.time()[ "elapsed" ]
  cat( "Beginning of reading data...\n" )

  # Read the two files: One contains the coordinates, and the other one 
  # contains the variables.
  strFileNameCoor = paste( "Data/", strTestName, "_coor.csv", sep = "" )
  coor            = read.csv( file = strFileNameCoor, sep = ";", head = F )
  strFileNameVars = paste( "Data/", strTestName, "_vars.csv", sep = "" )
  vars            = read.csv( file = strFileNameVars, sep = ";" )

  # Remove the first column of the variables since it is the line number.
  allVariables = vars[ , c( 2, 3, 4 ) ];
  rm( vars )

  # Set variables m, n, and k. 
  # Then, extract and organize coordinates into allCoordinates.
  m              = 3
  n              = dim( coor )[ 1 ]
  landmarks      = as.matrix( coor[ , 2 : ( dim( coor )[ 2 ] ) ] )
  k              = ( dim( landmarks )[ 2 ] ) / 3
  auxCoordinates = array( c( landmarks ), dim = c( n, 3, k ) )
  allCoordinates = aperm( auxCoordinates, c( 3, 2, 1 ) )
  rm( coor )
  rm( landmarks )
  rm( auxCoordinates )

  # Compute the object sizes.
  allSizes = apply( allCoordinates[ , , ], 3, my_new_centroid.size )
  #### write.table(sizes,'Data/sizes.txt')

  # Compute the preshapes of all objects.
  cat( "Beginning computation of preshapes...\n" )
  allPreshapes = array( rep( 0, n * ( m * k - m ) ), dim = c( k-1, m, n ) )
  for( i in seq( 1, n ) ) {
    allPreshapes[ , , i ] = my_new_preshape( allCoordinates[ , , i ] )
  }

  # Create the dataset: a list with the obtained info.
  dataset = list( k              = k,
                  n              = n,
                  m              = m,
                  allCoordinates = allCoordinates,
                  allVariables   = allVariables,
                  allSizes       = allSizes,
                  allPreshapes   = allPreshapes )

  cat( "Dataset is loaded: \n" )
  print_dataset_dimensions( dataset )

  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'read_dataset': ", t2 - t1, "\n\n" )

  # Return result.
  return( dataset )
}

