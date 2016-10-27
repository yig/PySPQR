# Python wrapper for SuiteSparseQR

This module wraps the [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html)
decomposition function for use with [SciPy](http://www.scipy.org).
This is Matlab's sparse `[Q,R,E] = qr()`.
Until now, no one wrapped that function of SuiteSparse.

# Usage

    M = scipy.sparse.rand( 10, 10, density = 0.1 )
    Q, R, E, rank = qr( M )
    ## Should be zero:
    print( abs( Q*R - M*permutation_from_E(E) ).sum() )

# Installation

Copy these two files next to your source code.

# Dependencies

* (SciPy/NumPy)[http://www.scipy.org]
* [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html) (`brew install suitesparse`)
* [cffi](http://cffi.readthedocs.io/) (`pip install cffi`)

# License

Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
