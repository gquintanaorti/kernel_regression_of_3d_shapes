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
# Module:  mdl_my_new_preshape.R
# =============================================================================
#
# This module contains several modifications of some functions provided by the
# noteworthy "shapes" package by Ian L. Dryden. These new implementations 
# have been accelerated and are much faster for medium and large datasets 
# than the original codes.
#
# This new code includes some commented lines used to time the new code and 
# to compare it with the original one. These lines are marked with 
# four "#" chars.
#
# =============================================================================


# =============================================================================
my_new_preshape = function( x ) {

    #### # Set timer.
    #### t1 = proc.time()[ "elapsed" ]

    if (is.complex(x)) {
        k <- nrow(as.matrix(x))
        h <- my_new_defh(k - 1)
        zstar <- x
        ztem <- h %*% zstar
        size <- sqrt(diag(Re(st(ztem) %*% ztem)))
        if (is.vector(zstar)) 
            z <- ztem/size
        if (is.matrix(zstar)) 
            z <- ztem %*% diag(1/size)
    }
    else {
        if (length(dim(x)) == 3) {
            k <- dim(x)[1]
            h <- my_new_defh(k - 1)
            n <- dim(x)[3]
            m <- dim(x)[2]
            z <- array(0, c(k - 1, m, n))
            for (i in 1:n) {
                z[, , i] <- h %*% x[, , i]
                size <- my_new_centroid.size(x[, , i])
                z[, , i] <- z[, , i]/size
            }
        }
        else {
            k <- nrow(as.matrix(x))

            #### # Set timer.
            #### t01 = proc.time()[ "elapsed" ]
            #### 
            #### h <- my_new_defh(k - 1)
            #### ztem <- h %*% x
            #### 
            #### # Print elapsed time.
            #### t02 = proc.time()[ "elapsed" ]
            #### cat( "  1. Elapsed time in hx: ", t02-t01, "\n" )
            #### 
            #### # Set timer.
            #### t01 = proc.time()[ "elapsed" ]

            # Use of my_new_compute_hx_implicitely to accelerate the code.
            ztem <- my_new_compute_hx_implicitely( x )

            #### # Print elapsed time.
            #### t02 = proc.time()[ "elapsed" ]
            #### cat( "  2. Elapsed time in hx: ", t02-t01, "\n" )
            #### 
            #### cat( "abs norm de diff: ", 
            ####      norm( ztem - ztem2, type = "F" ), "\n" )
            #### cat( "rel norm de diff: ", 
            ####      norm( ztem - ztem2, type = "F" ) / 
            ####      norm( ztem, type = "F" ), "\n" )
            #### cat( "  dim de h:    ", dim( h ), "\n" )
            #### cat( "  dim de x:    ", dim( x ), "\n" )
            #### cat( "  dim de ztem: ", dim( ztem ), "\n" )

            size <- my_new_centroid.size(x)
            z <- ztem/size
        }
    }
    #### # Print elapsed time.
    #### t2 = proc.time()[ "elapsed" ]
    #### cat( "Elapsed time in my_new_preshape (seconds): ", t2 - t1, "\n" )

    return( z )
}



# =============================================================================
my_new_centroid.size = function( x ) {
    if ((is.vector(x) == FALSE) && is.complex(x)) {
        k <- nrow(x)
        n <- ncol(x)
        tem <- array(0, c(k, 2, n))
        tem[, 1, ] <- Re(x)
        tem[, 2, ] <- Im(x)
        x <- tem
    }
    {
        if (length(dim(x)) == 3) {
            n <- dim(x)[3]
            sz <- rep(0, times = n)
            k <- dim(x)[1]
            h <- my_new_defh(k - 1)
            for (i in 1:n) {
                xh <- h %*% x[, , i]
                sz[i] <- sqrt(sum(diag(t(xh) %*% xh)))
            }
            sz
        }
        else {
            if (is.vector(x) && is.complex(x)) {
                x <- cbind(Re(x), Im(x))
            }
            k <- nrow(x)

            #### h <- my_new_defh(k - 1)
            #### xh <- h %*% x
            # Use of my_new_compute_hx_implicitely to accelerate the code.
            xh <- my_new_compute_hx_implicitely( x )

            size <- sqrt(sum(diag(t(xh) %*% xh)))
            size
        }
    }
}

# =============================================================================
my_new_defh = function( nrow ) {

    #### # Set timer.
    #### t1 = proc.time()[ "elapsed" ]

    k <- nrow
    h <- matrix( 0, k, k + 1 )
    if( nrow > 0 ) {
      for( j in seq( 1, k ) ) {
          val = -1 / sqrt( j * ( j + 1 ) )
          h[ j, seq( 1, j ) ] = val
          h[ j, j+1 ]         = - j * val
      }
    }

    #### # Print elapsed time.
    #### t2 = proc.time()[ "elapsed" ]
    #### cat( "Elapsed time in my_new_defh (seconds): ", t2 - t1, "\n" )

    return( h )
}

# =============================================================================
my_new_compute_hx_implicitely = function( x ) {

    #### # Set timer.
    #### t1 = proc.time()[ "elapsed" ]

    m = dim( x )[ 1 ]-1
    n = dim( x )[ 2 ]

    #### print( dim( x ) )

    result <- matrix( 0, m, n )
    vsum   <- rep( 0, n )
    if( m > 0 ) {
      for( i in seq( 1, m ) ) {
        val  = -1 / sqrt( i * ( i + 1 ) )
        hi   = val
        hip1 = - i * val
        #### cat( "i: ", i, "  hi:", hi, "  hip1: ", hip1, "\n" )

        #### for( j in seq( 1, n ) ) {
        ####   vsum[ j ]      = vsum[ j ] + x[ i, j ]
        ####   #### cat( "  j: ", j, " vsum[ j ]: ", vsum[ j ], "\n" )
        ####   result[ i, j ] = vsum[ j ] *  hi + x[ i+1, j ] * hip1
        #### }
        vsum = vsum + x[ i, ]
        result[ i, ] = vsum * hi + x[ i + 1, ] * hip1
      }
    }

    #### # Print elapsed time.
    #### t2 = proc.time()[ "elapsed" ]
    #### cat( "Elapsed time in my_new_defh (seconds): ", t2 - t1, "\n" )

    return( result )
}

