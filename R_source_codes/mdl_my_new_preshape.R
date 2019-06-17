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
my_preshape = function( x ) {
#
# Written by G. Quintana-Ortí and Amelia Simó, University Jaume I, Spain.
# It is an optimization for performance of Dryden's original code.
# This code is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED.
#
    if (is.complex(x)) {
        k <- nrow(as.matrix(x))
        h <- my_defh(k - 1)
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
            h <- my_defh(k - 1)
            n <- dim(x)[3]
            m <- dim(x)[2]
            z <- array(0, c(k - 1, m, n))
            for (i in 1:n) {
                z[, , i] <- h %*% x[, , i]
                size <- my_centroid.size(x[, , i])
                z[, , i] <- z[, , i]/size
            }
        }
        else {
            k <- nrow(as.matrix(x))

            #### h <- my_defh(k - 1)
            #### ztem <- h %*% x

            # Use of multiply_by_helmert_implicitly to accelerate the code.
            ztem <- multiply_by_helmert_implicitly( x )

            size <- my_centroid.size(x)
            z <- ztem/size
        }
    }

    z
}

# =============================================================================
my_centroid.size = function( x ) {
#
# Written by G. Quintana-Ortí and Amelia Simó, University Jaume I, Spain.
# It is an optimization for performance of Dryden's original code.
# This code is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED.
#
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
            h <- my_defh(k - 1)
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

            #### h <- my_defh(k - 1)
            #### xh <- h %*% x

            # Use of multiply_by_helmert_implicitly to accelerate the code.
            xh <- multiply_by_helmert_implicitly( x )

            size <- sqrt(sum(diag(t(xh) %*% xh)))
            size
        }
    }
}

# =============================================================================
my_defh = function( nrow ) {
#
# Written by G. Quintana-Ortí and Amelia Simó, University Jaume I, Spain.
# It is an optimization for performance of Dryden's original code.
# This code is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED.
#
    k <- nrow
    h <- matrix( 0, k, k + 1 )
    if( nrow > 0 ) {
      for( j in seq( 1, k ) ) {
          val = -1 / sqrt( j * ( j + 1 ) )
          h[ j, seq( 1, j ) ] = val
          h[ j, j+1 ]         = - j * val
      }
    }

    h
}

# =============================================================================
multiply_by_helmert_implicitly = function( x ) {
#
# This code multiplies the "x" argument by the Helmert matrix of the 
# corresponding size.
# Written by G. Quintana-Ortí and Amelia Simó, University Jaume I, Spain.
# This code is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED.
#
    m = dim( x )[ 1 ]-1
    n = dim( x )[ 2 ]

    result <- matrix( 0, m, n )
    vsum   <- rep( 0, n )
    if( m > 0 ) {
      for( i in seq( 1, m ) ) {
        val  = -1 / sqrt( i * ( i + 1 ) )
        hi   = val
        hip1 = - i * val

        vsum = vsum + x[ i, ]
        result[ i, ] = vsum *  hi + x[ i+1, ] * hip1

      }
    }

    return( result )
}


