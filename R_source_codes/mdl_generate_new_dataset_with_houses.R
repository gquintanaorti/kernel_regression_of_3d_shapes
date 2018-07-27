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
# Module:  mdl_generate_new_dataset_with_houses.R
# =============================================================================
#
# This module contains several functions related to the generation of a new
# dataset with houses.
#
# =============================================================================


# =============================================================================
generate_files_with_dataset_of_houses = 
    function( numLandmarksPerEdge, errorFactor ) {
#
# It overwrites two files ("Data/houses1_coor.csv" and "Data/houses1_vars.csv")
# with new data on several (400) houses.
# The first output file contains the coordinates. The second output file 
# contains the set of variables (varset) employed to define the houses.
#
# Method's arguments:
#   numLandmarksPerEdge:  Number of landmarks in every edge.
#   errorFactor:          Factor applied to the landmarks of the houses.
#

  # Define the number of houses per every varset (set of variables).
  # Since there are 8 varsets, the total number of houses will 
  # be 8 * numHousesPerVarset.
  numHousesPerVarset = 50

  # Seet seed for random number generator. 
  # A fixed initial seed is set for reproducibility purposes. If different 
  # data wants to be obtained, next line should be commented out.
  set.seed( 1 )

  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # Define the two file names.
  strFileNameCoor = "Data/houses1_coor.csv"
  strFileNameVars = "Data/houses1_vars.csv"

  # Remove the files before appending data into them.
  unlink( strFileNameCoor )
  unlink( strFileNameVars )

  # Generate header for file with vars.
  # File with coordinates (landmarks) do not include a header since the number
  # of landmarks could be very large.
  newvarsetRow = data.frame( 
                     as.list( c( "lineNumber", "varX", "varY", "varZ" ) ) )
  save_data_into_file( strFileNameVars, newvarsetRow )

  # Generate the initial model house.
  modelHouse = generate_model_house( numLandmarksPerEdge )

  # Generate several types of houses (with different sets of variables).
  lineNumber = 1
  for( x in c( 10, 20 ) ) {
    for( y in c( 10, 20 ) ) {
      for( z in c( 10, 20 ) ) {

        # Generate multiple houses for every different varset.
        for( i in seq( 1, numHousesPerVarset ) ) {
          # Generate and save one house into both files.
          varset       = c( x, y, z )
          house        = generate_a_house( modelHouse, varset, errorFactor )
          newCoordRow  = data.frame( 
                             as.list( c( lineNumber, c( t( house ) ) ) ) )
          newvarsetRow = data.frame( 
                             as.list( c( lineNumber, c( varset ) ) ) )
          save_data_into_file( strFileNameCoor, newCoordRow )
          save_data_into_file( strFileNameVars, newvarsetRow )
          lineNumber = lineNumber + 1
        }
      }
    }
  }

  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'generate_files_with_dataset_of_houses': ", 
       t2 - t1, "\n\n" )
}

# =============================================================================
generate_model_house = function( numLandmarksPerEdge ) {
#
# It builds an initial model house with "n" vertices along each edge.
#
  # Check the number of landmarks along each edge.
  if( numLandmarksPerEdge < 2 ) {
    cat( "ERROR in generate_model_house: numLandmarksPerEdge is too small.\n" )
    cat( "numLandmarksPerEdge should be >= 2.\n" )
    stop()
  }

  # Build and save the vertices of the initial house into a matrix.
  verticesOfModelHouse = matrix( c(   0,   0,   0,
                                      1,   0,   0,
                                      1, 0.7,   0,
                                    0.5, 1.0,   0,
                                      0, 0.7,   0,
                                      0,   0,   1,
                                      1,   0,   1,
                                      1, 0.7,   1,
                                    0.5, 1.0,   1,
                                      0, 0.7,   1 ), byrow = TRUE, ncol = 3 )

  # Build the edges of the initial house by generating multiple intermediate
  # landmarks between some vertices.
  modelHouse = verticesOfModelHouse
  if( numLandmarksPerEdge > 2 ) {
    modelHouse = rbind( modelHouse,
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  1,  2 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  2,  3 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  3,  4 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  4,  5 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  5,  1 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  6,  7 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  7,  8 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  8,  9 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  9, 10 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse, 10,  6 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  1,  6 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  2,  7 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  3,  8 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  4,  9 ),
        generate_edge( numLandmarksPerEdge, verticesOfModelHouse,  5, 10 ) )
  }

  # Return results.
  #### print( modelHouse )
  return( modelHouse )
}

# =============================================================================
generate_a_house = function( modelHouse, varset, errorFactor ) {
#
# It generates a new house based on "modelHouse" by applying a scaling defined
# by varset to every coordinate. 
# "varset" must be a vector with 3 components (variables).
# It returns the matrix with the coordinates.
#
  # Create new house.
  newHouse = modelHouse

  # Add some error to all landmarks.
  numElems = nrow( modelHouse ) * ncol( modelHouse )
  newHouse = newHouse + rnorm( n = numElems, sd = errorFactor )

  # Scale new house according to varset.
  newHouse = newHouse %*% diag( varset )

  # Round coordinates of new house.
  newHouse = round( newHouse, 2 )

  # Return results.
  return( newHouse )
}

# =============================================================================
display_house = function( house ) {
#
# It displays the house on the screen. All landmarks are shown as points, but 
# several lines are added to better show the shape of the house. 
# It saves the image in the "house1.pdf" file.
#
  # Libraries to be used.
  library( scatterplot3d )
  library( rgl )

  # Shift house to origin.
  house = shift_to_origin( house )

  # Plot the landmarks with points.
  plot3d( house, 
          aspect = F, axes = F,
          xlab = "", ylab = "", zlab = "", 
          size = 10, type = "p" )

  # Plot the edges with lines.
  plot3d( house[ c(  1,  2 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  2,  3 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  3,  4 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  4,  5 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  5,  1 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  6,  7 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  7,  8 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  8,  9 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  9, 10 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c( 10,  6 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  1,  6 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  2,  7 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  3,  8 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  4,  9 ), ], aspect = F, type = "l", add = T )
  plot3d( house[ c(  5, 10 ), ], aspect = F, type = "l", add = T )

  # Save image into a file.
  rgl.postscript( "house1.pdf", "pdf" )
}

# =============================================================================
print_images_of_houses_without_error = function() {
#
# It prints the images of the 8 houses without error into EPS files.
# This method can be used to generate and print the eight types of houses.
# Employed to generate images of the houses for the paper.
#
  # Libraries to be used.
  library( scatterplot3d )
  library( rgl )

  # Set timer.
  t1 = proc.time()[ "elapsed" ]

  # Generate the initial house.
  modelHouse = generate_model_house( 3 )
  rgl.postscript( "original_house.eps" )
  rgl.postscript( "original_house.pdf", "pdf" )

  # Generate and save one house type.
  for( x in c( 10, 20 ) ) {
    for( y in c( 10, 20 ) ) {
      for( z in c( 10, 20 ) ) {
        varset = c( x, y, z )
        house  = generate_a_house( modelHouse, varset, 0.0 )
        display_house( house )
        cat( "\n" )

        # Save in EPS format.
        fileName = paste( "house_x", x, "_y", y, "_z", z, ".eps", sep = "" )
        rgl.postscript( fileName )

        # Save in PDF format.
        fileName = paste( "house_x", x, "_y", y, "_z", z, ".pdf", sep = "" )
        rgl.postscript( fileName, "pdf" )
      }
    }
  }

  # Print elapsed time.
  t2 = proc.time()[ "elapsed" ]
  cat( "Elapsed time (s) in 'print_images_of_houses_without_error': ", 
       t2 - t1, "\n\n" )
}

# =============================================================================
save_data_into_file = function( strFileName, rowWithData ) {
#
# It adds data received in arguments into file.
#
  write.table( format( rowWithData, decimal.mark = "." ), file = strFileName,
      sep = ";", quote = FALSE, row.names = FALSE, col.names = FALSE,
      append = TRUE )
}

# =============================================================================
generate_edge = function( numLandmarksPerEdge, house, idx1, idx2 ) {
#
# It generates a discrete edge, that is, a sequence of numLandmarksPerEdge-2
# landmarks between "idx1" landmark and "idx2" landmark of "house".
#
  edge = generate_intermediate_landmarks( numLandmarksPerEdge, 
             as.vector( house[ idx1, ] ),
             as.vector( house[ idx2, ] ) )
  # Return results.
  return( edge )
}

# =============================================================================
generate_intermediate_landmarks = function( numLandmarksPerEdge, a, b ) {
#
# It generates a matrix with numLandmarksPerEdge-2 intermediate landmarks 
# between landmarks "a" and "b", without including both ends.
#
  # Check input arguments.
  if( numLandmarksPerEdge < 2 ) {
    cat( "ERROR in generate_intermediate_landmarks: " )
    cat( "Invalid numLandmarksPerEdge: ", numLandmarksPerEdge, "\n" )
    stop()
  }

  # Generate intermediate landmarks.
  vx = seq( a[ 1 ], b[ 1 ], length.out = numLandmarksPerEdge )
  vy = seq( a[ 2 ], b[ 2 ], length.out = numLandmarksPerEdge )
  vz = seq( a[ 3 ], b[ 3 ], length.out = numLandmarksPerEdge )
  vecLandmarks = c( vx, vy, vz )
  landmarks = matrix( vecLandmarks, byrow = F, ncol = 3 )

  # Remove initial and final landmarks.
  landmarks = landmarks[ seq( 2, numLandmarksPerEdge - 1 ), ]

  # Return results.
  return( landmarks )
}

# =============================================================================
center_matrix = function( a ) {
#
# It centers the contents of data for every column.
#

  b     = a
  ncols = dim( a )[ 2 ]
  for( j in seq( 1, ncols ) ) {
    colMean = mean( a[ , j ] )
    b[ , j ] = a[ , j ] - colMean
  }

  # Return results.
  return( b )
}

# =============================================================================
shift_to_origin = function( a ) {
#
# It shifts the matrix to the origin by substracting the minimun of every
# column to every column.
#

  b     = a
  ncols = dim( a )[ 2 ]
  for( j in seq( 1, ncols ) ) {
    colMin = min( a[ , j ] )
    b[ , j ] = a[ , j ] - colMin
    cat( "  Maximum value in col: ", j, " is: ", max( b[ , j ] ), "\n" )
  }

  # Return results.
  return( b )
}

