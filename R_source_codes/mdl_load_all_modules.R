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
# Module:  mdl_load_all_modules.R
# =============================================================================
#
# This module removes all the code and then loads the code contained in all 
# the rest of modules.
#
# =============================================================================

#
# Set the "locale" variable to remove some warnings on Linux systems.
#
Sys.setlocale( 'LC_ALL', 'C' )

#
# Define "ls_fun" function, just in case it is not defined.
#
# =============================================================================
ls_fun <- function( name = parent.frame() ) {
  obj <- ls( name = name )
  obj[ sapply( obj, function( x ) is.function( get( x ) ) ) ]
}

#
# Remove all the previously loaded code, but keep all the data.
#
rm( list = ls_fun() )

#
# Redefine "ls_no_fun" and "ls_fun" functions because all the code has been
# removed.
#
# =============================================================================
ls_no_fun <- function( name = parent.frame() ) {
  obj <- ls( name = name )
  obj[ ! sapply( obj, function( x ) is.function( get( x ) ) ) ]
}

# =============================================================================
ls_fun <- function( name = parent.frame() ) {
  obj <- ls( name = name )
  obj[ sapply( obj, function( x ) is.function( get( x ) ) ) ]
}

#
# Load all the modules with R code.
#
source( "mdl_generate_new_dataset_with_houses.R" )
source( "mdl_read_dataset.R" )
source( "mdl_dataset_utils.R" )
source( "mdl_my_new_preshape.R" )
source( "mdl_kernel_regression_in_shape_space.R" )
source( "mdl_tests.R" )
source( "mdl_distances.R" )
source( "mdl_cross_validation.R" )

