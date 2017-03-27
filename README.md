# Python wrapper for SuiteSparseQR

This module wraps the [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html)
decomposition function for use with [SciPy](http://www.scipy.org).
This is Matlab's sparse `[Q,R,E] = qr()`.
For some reason, no one ever wrapped that function of SuiteSparseQR for Python.

# Usage

    # Perform QR decomposition of a sparse matrix
    import spqr
    import scipy.sparse
    M = scipy.sparse.rand( 10, 10, density = 0.1 )
    Q, R, E, rank = spqr.qr( M )
    ## Should be approximately zero:
    print( abs( Q*R - M*spqr.permutation_from_E(E) ).sum() )

    # Solve an overdetermined linear system  A x = b  in the least-squares sense
    #
    # (The same routine also works for the usual non-overdetermined case.)
    #
    A = scipy.sparse.rand( 20, 10, density = 0.1 )  # 20 equations, 10 unknowns
    b = numpy.random.random(10)  # one RHS, but could also have many (in shape (10,k))
    x = qr_solve( A, b, tolerance = 0 )

    # Solve many linear systems "M x = b for b in columns(B)"
    #
    B = scipy.sparse.rand( 10, 5, density = 0.1 )  # many RHS, sparse
    x = qr_solve_sparse( M, B, tolerance = 0 )

    # Solve a system  M x = b  via QR decomposition
    #
    # This approach is slow due to explicit construction of Q, but may be
    # useful if a large number of systems need to be solved with the same M.
    #
    import numpy
    import scipy
    Q, R, E, rank = spqr.qr( M )
    k = min(M.shape)
    R = R.tocsr()[:k,:]
    Q = Q.tocsc()[:,:k]
    Qb = (Q.T).dot(b)
    result = scipy.sparse.linalg.spsolve(R, Qb)
    x = numpy.empty_like(result)
    x[E] = result[:]


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
