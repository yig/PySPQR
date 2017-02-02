# Python wrapper for SuiteSparseQR

This module wraps the [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html)
decomposition function for use with [SciPy](http://www.scipy.org).
This is Matlab's sparse `[Q,R,E] = qr()`.
For some reason, no one ever wrapped that function of SuiteSparseQR for Python.

# Usage

    import spqr
    import scipy.sparse
    M = scipy.sparse.rand( 10, 10, density = 0.1 )
    Q, R, E, rank = spqr.qr( M )
    ## Should be approximately zero:
    print( abs( Q*R - M*spqr.permutation_from_E(E) ).sum() )

# Installation

Copy the three `.py` files next to your source code,
or leave them in this directory and call it as a module.

Tested on Python 2.7 and 3.5.

Tested on Mac OS X and Ubuntu Linux.

# Dependencies

* [SciPy/NumPy](http://www.scipy.org)
* [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html) (`brew install suitesparse`)
* [cffi](http://cffi.readthedocs.io/) (`pip install cffi`)

# License

Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
