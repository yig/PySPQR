'''
Author: Yotam Gingold <yotam (strudel) yotamgingold.com>
License: Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
Description: Wrapper for SuiteSparse qr() function. Matlab has it, Python should have it, too.
'''

from __future__ import print_function, division

try:
    from _spqr import ffi, lib
except ImportError:
    print( "=== Wrapper module not compiled; compiling..." )
    import spqr_gen
    spqr_gen.main()
    print( "=== ...compiled." )

    from _spqr import ffi, lib

import scipy.sparse
import numpy
import cffi_asarray

'''
Helpful links:
    The primary docs:
    http://cffi.readthedocs.io/en/latest/overview.html
    http://cffi.readthedocs.io/en/latest/using.html

    Some helpful examples for #define and NumPy:
    https://kogs-www.informatik.uni-hamburg.de/~seppke/content/teaching/wise1314/20131107_pridoehl-cffi.pdf
'''

## Initialize cholmod
cc = ffi.new("cholmod_common*")
lib.cholmod_l_start( cc )

## Set up cholmod deinit to run when Python exits
def _deinit():
    '''Deinitialize the CHOLMOD library.'''
    lib.cholmod_l_finish( cc )
import atexit
atexit.register(_deinit)


## Data format conversion

def scipy2cholmodsparse( scipy_A ):
    '''Convert a SciPy sparse matrix to a CHOLMOD sparse matrix.

The input is first internally converted to scipy.sparse.coo_matrix format.

When no longer needed, the returned CHOLMOD sparse matrix must be deallocated using cholmod_free_sparse().
'''
    scipy_A = scipy_A.tocoo()

    nnz = scipy_A.nnz

    ## There is a potential performance win if we know A is symmetric and we
    ## can get only the upper or lower triangular elements.
    chol_A = lib.cholmod_l_allocate_triplet( scipy_A.shape[0], scipy_A.shape[1], nnz, 0, lib.CHOLMOD_REAL, cc )

    Ai = ffi.cast( "SuiteSparse_long*", chol_A.i )
    Aj = ffi.cast( "SuiteSparse_long*", chol_A.j )
    Avals = ffi.cast( "double*", chol_A.x )

    Ai[0:nnz] = scipy_A.row
    Aj[0:nnz] = scipy_A.col
    Avals[0:nnz] = scipy_A.data

    chol_A.nnz = nnz

    ## Print what cholmod sees as the matrix.
    # lib.cholmod_l_print_triplet( chol_A, "A".encode('utf-8'), cc )
    assert lib.cholmod_l_check_triplet( chol_A, cc ) == 1

    ## Convert to a cholmod_sparse matrix.
    result = lib.cholmod_l_triplet_to_sparse( chol_A, nnz, cc )
    ## Free the space used by the cholmod triplet matrix.
    _cholmod_free_triplet( chol_A )

    return result

def cholmodsparse2scipy( chol_A ):
    '''Convert a CHOLMOD sparse matrix to a scipy.sparse.coo_matrix.'''
    ## Convert to a cholmod_triplet matrix.
    chol_A = lib.cholmod_l_sparse_to_triplet( chol_A, cc )

    nnz = chol_A.nnz

    Ai = ffi.cast( "SuiteSparse_long*", chol_A.i )
    Aj = ffi.cast( "SuiteSparse_long*", chol_A.j )
    Adata = ffi.cast( "double*", chol_A.x )

    ## Have to pass through list().
    ## https://bitbucket.org/cffi/cffi/issues/292/cant-copy-data-to-a-numpy-array
    ## http://stackoverflow.com/questions/16276268/how-to-pass-a-numpy-array-into-a-cffi-function-and-how-to-get-one-back-out
    '''
    i = numpy.zeros( nnz, dtype = numpy.int64 )
    j = numpy.zeros( nnz, dtype = numpy.int64 )
    data = numpy.zeros( nnz )

    i[0:nnz] = list( Ai[0:nnz] )
    j[0:nnz] = list( Aj[0:nnz] )
    data[0:nnz] = list( Adata[0:nnz] )
    '''
    ## UPDATE: I can do this without going through list() or making two extra copies.
    ## NOTE: Create a copy() of the array data, because the coo_matrix() constructor
    ##       doesn't and the cholmod memory fill get freed.
    i = cffi_asarray.asarray( ffi, Ai, nnz ).copy()
    j = cffi_asarray.asarray( ffi, Aj, nnz ).copy()
    data = cffi_asarray.asarray( ffi, Adata, nnz ).copy()

    scipy_A = scipy.sparse.coo_matrix(
        ( data, ( i, j ) ),
        shape = ( chol_A.nrow, chol_A.ncol )
        )

    ## Free the space used by the cholmod triplet matrix.
    _cholmod_free_triplet( chol_A )

    return scipy_A

def numpy2cholmoddense( numpy_A ):
    '''Convert a NumPy array (rank-1 or rank-2) to a CHOLMOD dense matrix.

Rank-1 arrays are converted to column vectors.

When no longer needed, the returned CHOLMOD sparse matrix must be deallocated using cholmod_free_dense().
'''
    numpy_A = numpy.atleast_2d( numpy_A )
    if numpy_A.shape[0] == 1 and numpy_A.shape[1] > 1:  # prefer column vector
        numpy_A = numpy_A.T
    nrow = numpy_A.shape[0]
    ncol = numpy_A.shape[1]
    lda  = nrow  # cholmod_dense is column-oriented
    chol_A = lib.cholmod_l_allocate_dense( nrow, ncol, lda, lib.CHOLMOD_REAL, cc )
    if chol_A == ffi.NULL:
        raise RuntimeError("Failed to allocate chol_A")
    Adata = ffi.cast( "double*", chol_A.x )
    for j in range(ncol):  # FIXME inefficient?
        Adata[(j*lda):((j+1)*lda)] = numpy_A[:,j]
    return chol_A

def cholmoddense2numpy( chol_A ):
    '''Convert a CHOLMOD dense matrix to a NumPy array.'''
    Adata = ffi.cast( "double*", chol_A.x )

    result = cffi_asarray.asarray( ffi, Adata, chol_A.nrow*chol_A.ncol )
    result = result.reshape( (chol_A.nrow, chol_A.ncol), order='F' )
    return numpy.squeeze(result)

def permutation_from_E( E ):
    '''Convert a permutation vector E (list or rank-1 array, length n) to a permutation matrix (n by n).

The result is returned as a scipy.sparse.coo_matrix, where the entries at (E[k], k) are 1.
'''
    n = len( E )
    j = numpy.arange( n )
    return scipy.sparse.coo_matrix( ( numpy.ones(n), ( E, j ) ), shape = ( n, n ) )


## Memory management

# Used only internally by this module (the user sees only sparse and dense formats).
def _cholmod_free_triplet( A ):
    '''Deallocate a CHOLMOD triplet format matrix.'''
    A_ptr = ffi.new("cholmod_triplet**")
    A_ptr[0] = A
    lib.cholmod_l_free_triplet( A_ptr, cc )

def cholmod_free_sparse( A ):
    '''Deallocate a CHOLMOD sparse matrix.'''
    A_ptr = ffi.new("cholmod_sparse**")
    A_ptr[0] = A
    lib.cholmod_l_free_sparse( A_ptr, cc )

def cholmod_free_dense( A ):
    '''Deallocate a CHOLMOD dense matrix.'''
    A_ptr = ffi.new("cholmod_dense**")
    A_ptr[0] = A
    lib.cholmod_l_free_dense( A_ptr, cc )


## Solvers

def qr( A, tolerance = None ):
    '''
    Given a sparse matrix A,
    returns Q, R, E, rank such that:
        Q*R = A*permutation_from_E(E)
    rank is the estimated rank of A.

    If optional `tolerance` parameter is negative, it has the following meanings:
        #define SPQR_DEFAULT_TOL ...       /* if tol <= -2, the default tol is used */
        #define SPQR_NO_TOL ...            /* if -2 < tol < 0, then no tol is used */

    The performance-optimal format for A is scipy.sparse.coo_matrix.

    For solving systems of the form A x = b, see qr_solve().

    qr() can also be used to solve systems, as follows:

        # inputs: scipy.sparse.coo_matrix A, rank-1 numpy.array b (RHS)
        import numpy
        import scipy

        Q, R, E, rank = spqr.qr( A )
        k = min(A.shape)
        R = R.tocsr()[:k,:]  # cut away the block of zeros at the bottom
        Q = Q.tocsc()[:,:k]  # cut away the corresponding part of Q
        Qb = (Q.T).dot(b)    # use .dot() (not numpy.dot()) since Q is sparse
        result = scipy.sparse.linalg.spsolve(R, Qb)
        x = numpy.empty_like(result)
        x[E] = result[:]     # apply inverse permutation

    Be aware that this approach is slow, because qr() explicitly constructs Q.
    Unless you have a large number of systems to solve with the same A, qr_solve() is much faster.
    '''

    chol_A = scipy2cholmodsparse( A )

    chol_Q = ffi.new("cholmod_sparse**")
    chol_R = ffi.new("cholmod_sparse**")
    chol_E = ffi.new("SuiteSparse_long**")

    if tolerance is None: tolerance = lib.SPQR_DEFAULT_TOL

    rank = lib.SuiteSparseQR_C_QR(
        ## Input
        lib.SPQR_ORDERING_DEFAULT,
        tolerance,
        A.shape[0],
        chol_A,
        ## Output
        chol_Q,
        chol_R,
        chol_E,
        cc
        )

    scipy_Q = cholmodsparse2scipy( chol_Q[0] )
    scipy_R = cholmodsparse2scipy( chol_R[0] )

    ## If chol_E is null, there was no permutation.
    if chol_E == ffi.NULL:
        E = None
    else:
        ## Have to pass through list().
        ## https://bitbucket.org/cffi/cffi/issues/292/cant-copy-data-to-a-numpy-array
        ## http://stackoverflow.com/questions/16276268/how-to-pass-a-numpy-array-into-a-cffi-function-and-how-to-get-one-back-out
        # E = numpy.zeros( A.shape[1], dtype = int )
        # E[0:A.shape[1]] = list( chol_E[0][0:A.shape[1]] )
        ## UPDATE: I can do this without going through list() or making two extra copies.
        E = cffi_asarray.asarray( ffi, chol_E[0], A.shape[1] ).copy()

    ## Free cholmod stuff
    cholmod_free_sparse( chol_Q[0] )
    cholmod_free_sparse( chol_R[0] )
    ## Apparently we don't need to do this. (I get a malloc error.)
    # lib.cholmod_l_free( A.shape[1], ffi.sizeof("SuiteSparse_long"), chol_E, cc )

    return scipy_Q, scipy_R, E, rank

def qr_solve( A, b, tolerance = None ):
    '''
    Given a sparse m-by-n matrix A, and dense m-by-k matrix (storing k RHS vectors) b,
    solve A x = b in the least-squares sense.

    This is much faster than using qr() to solve the system, since Q is not explicitly constructed.

    Returns x on success, None on failure.

    The performance-optimal format for A is scipy.sparse.coo_matrix.

    If optional `tolerance` parameter is negative, it has the following meanings:
        #define SPQR_DEFAULT_TOL ...       /* if tol <= -2, the default tol is used */
        #define SPQR_NO_TOL ...            /* if -2 < tol < 0, then no tol is used */
    '''

    chol_A = scipy2cholmodsparse( A )
    chol_b = numpy2cholmoddense(  b )

    if tolerance is None: tolerance = lib.SPQR_DEFAULT_TOL

    chol_x = lib.SuiteSparseQR_C_backslash(
        ## Input
        lib.SPQR_ORDERING_DEFAULT,
        tolerance,
        chol_A,
        chol_b,
        cc )

    if chol_x == ffi.NULL:
        return None  # failed

    numpy_x = cholmoddense2numpy( chol_x )

    ## Free cholmod stuff
    cholmod_free_sparse( chol_A )
    cholmod_free_dense(  chol_b )
    cholmod_free_dense(  chol_x )

    return numpy_x


if __name__ == '__main__':
    # Q, R, E, rank = qr( scipy.sparse.identity(10) )

    M = scipy.sparse.rand( 10, 10, density = 0.1 )
    Q, R, E, rank = qr( M, tolerance = 0 )
    print( abs( Q*R - M*permutation_from_E(E) ).sum() )

    b = numpy.random.random(10)  # one RHS, but could also have many (in shape (10,k))
    x = qr_solve( M, b, tolerance = 0 )
    print( x )
