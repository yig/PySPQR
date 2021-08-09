'''
Author: Yotam Gingold <yotam (strudel) yotamgingold.com>
License: Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
Description: Wrapper for SuiteSparse qr() and solve() functions. Matlab and Julia have it, Python should have it, too.
'''

from __future__ import print_function, division, absolute_import
import os
from os.path import join, expanduser
import platform

from cffi import FFI

include_dirs = [ '/usr/include/suitesparse', join('C:', 'Program Files', 'Python', 'suitesparse') ]
libraries = ['spqr']

# for compatibility with conda envs
if 'CONDA_DEFAULT_ENV' in os.environ:
    homedir = expanduser("~")
    include_dirs.append( join(homedir, 'anaconda3', 'envs', os.environ['CONDA_DEFAULT_ENV'], 'Library', 'include', 'suitesparse') )
    include_dirs.append( join(homedir, 'miniconda3', 'envs', os.environ['CONDA_DEFAULT_ENV'], 'Library', 'include', 'suitesparse') )

if platform.system() == 'Windows':
    # https://github.com/yig/PySPQR/issues/6
    libraries.extend( ['amd','btf','camd','ccolamd','cholmod','colamd','cxsparse'
'klu','lapack','ldl','lumfpack','metis','suitesparseconfig','libblas'] )

ffibuilder = FFI()

ffibuilder.set_source( "sparseqr._sparseqr",
    """#include <SuiteSparseQR_C.h>
""",
    ## You may need to modify the following line,
    ## which is needed on Ubuntu and harmless on Mac OS.
    include_dirs = include_dirs,
    libraries = libraries )

ffibuilder.cdef("""
// The int... is a magic thing which tells the compiler to figure out what the right
// integer type is.
typedef int... SuiteSparse_long;

/// Many of these are copied from "cholmod_core.h"

// The cholmod_common struct can't be completely opaque,
// since we need to allocate space for one.
typedef struct cholmod_common { ...; } cholmod_common ;

// We can keep the cholmod_sparse struct opaque since we will only ever
// interact with it by converting to and from a triplet struct.
typedef ... cholmod_sparse ;
typedef struct cholmod_triplet_struct
{
    size_t nrow ;	/* the matrix is nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;	/* maximum number of entries in the matrix */
    size_t nnz ;	/* number of nonzeros in the matrix */

    void *i ;		/* i [0..nzmax-1], the row indices */
    void *j ;		/* j [0..nzmax-1], the column indices */
    void *x ;		/* size nzmax or 2*nzmax, if present */
    void *z ;		/* size nzmax, if present */

    int stype ;		/* Describes what parts of the matrix are considered:
			 *
	* 0:  matrix is "unsymmetric": use both upper and lower triangular parts
	*     (the matrix may actually be symmetric in pattern and value, but
	*     both parts are explicitly stored and used).  May be square or
	*     rectangular.
	* >0: matrix is square and symmetric.  Entries in the lower triangular
	*     part are transposed and added to the upper triangular part when
	*     the matrix is converted to cholmod_sparse form.
	* <0: matrix is square and symmetric.  Entries in the upper triangular
	*     part are transposed and added to the lower triangular part when
	*     the matrix is converted to cholmod_sparse form.
	*
	* Note that stype>0 and stype<0 are different for cholmod_sparse and
	* cholmod_triplet.  The reason is simple.  You can permute a symmetric
	* triplet matrix by simply replacing a row and column index with their
	* new row and column indices, via an inverse permutation.  Suppose
	* P = L->Perm is your permutation, and Pinv is an array of size n.
	* Suppose a symmetric matrix A is represent by a triplet matrix T, with
	* entries only in the upper triangular part.  Then the following code:
	*
	*	Ti = T->i ;
	*	Tj = T->j ;
	*	for (k = 0 ; k < n  ; k++) Pinv [P [k]] = k ;
	*	for (k = 0 ; k < nz ; k++) Ti [k] = Pinv [Ti [k]] ;
	*	for (k = 0 ; k < nz ; k++) Tj [k] = Pinv [Tj [k]] ;
	*
	* creates the triplet form of C=P*A*P'.  However, if T initially
	* contains just the upper triangular entries (T->stype = 1), after
	* permutation it has entries in both the upper and lower triangular
	* parts.  These entries should be transposed when constructing the
	* cholmod_sparse form of A, which is what cholmod_triplet_to_sparse
	* does.  Thus:
	*
	*	C = cholmod_triplet_to_sparse (T, 0, &Common) ;
	*
	* will return the matrix C = P*A*P'.
	*
	* Since the triplet matrix T is so simple to generate, it's quite easy
	* to remove entries that you do not want, prior to converting T to the
	* cholmod_sparse form.  So if you include these entries in T, CHOLMOD
	* assumes that there must be a reason (such as the one above).  Thus,
	* no entry in a triplet matrix is ever ignored.
	*/

    int itype ; /* CHOLMOD_LONG: i and j are SuiteSparse_long.  Otherwise int */
    int xtype ; /* pattern, real, complex, or zomplex */
    int dtype ; /* x and z are double or float */

} cholmod_triplet ;


/* -------------------------------------------------------------------------- */
/* cholmod_start:  first call to CHOLMOD */
/* -------------------------------------------------------------------------- */
int cholmod_l_start (cholmod_common *) ;
/* -------------------------------------------------------------------------- */
/* cholmod_finish:  last call to CHOLMOD */
/* -------------------------------------------------------------------------- */
int cholmod_l_finish (cholmod_common *) ;

/* -------------------------------------------------------------------------- */
/* cholmod_free_sparse:  free a sparse matrix */
/* -------------------------------------------------------------------------- */
int cholmod_l_free_sparse
(
    /* ---- in/out --- */
    cholmod_sparse **A,	/* matrix to deallocate, NULL on output */
    /* --------------- */
    cholmod_common *Common
) ;

/* ========================================================================== */
/* === Core/cholmod_dense =================================================== */
/* ========================================================================== */

/* A dense matrix in column-oriented form.  It has no itype since it contains
 * no integers.  Entry in row i and column j is located in x [i+j*d].
 */

typedef struct cholmod_dense_struct
{
    size_t nrow ;	/* the matrix is nrow-by-ncol */
    size_t ncol ;
    size_t nzmax ;	/* maximum number of entries in the matrix */
    size_t d ;		/* leading dimension (d >= nrow must hold) */
    void *x ;		/* size nzmax or 2*nzmax, if present */
    void *z ;		/* size nzmax, if present */
    int xtype ;		/* pattern, real, complex, or zomplex */
    int dtype ;		/* x and z double or float */

} cholmod_dense ;

/* -------------------------------------------------------------------------- */
/* cholmod_allocate_dense:  allocate a dense matrix (contents uninitialized) */
/* -------------------------------------------------------------------------- */

cholmod_dense *cholmod_l_allocate_dense
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of matrix */
    size_t ncol,	/* # of columns of matrix */
    size_t d,		/* leading dimension */
    int xtype,		/* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_free_dense:  free a dense matrix */
/* -------------------------------------------------------------------------- */

int cholmod_l_free_dense
(
    /* ---- in/out --- */
    cholmod_dense **X,	/* dense matrix to deallocate, NULL on output */
    /* --------------- */
    cholmod_common *Common
) ;

/*
 * ============================================================================
 * === cholmod_triplet ========================================================
 * ============================================================================
 *
 * A sparse matrix held in triplet form is the simplest one for a user to
 * create.  It consists of a list of nz entries in arbitrary order, held in
 * three arrays: i, j, and x, each of length nk.  The kth entry is in row i[k],
 * column j[k], with value x[k].  There may be duplicate values; if A(i,j)
 * appears more than once, its value is the sum of the entries with those row
 * and column indices.
 */
/* -------------------------------------------------------------------------- */
/* cholmod_allocate_triplet:  allocate a triplet matrix */
/* -------------------------------------------------------------------------- */
cholmod_triplet *cholmod_l_allocate_triplet
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of T */
    size_t ncol,	/* # of columns of T */
    size_t nzmax,	/* max # of nonzeros of T */
    int stype,		/* stype of T */
    int xtype,		/* CHOLMOD_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
) ;
#define CHOLMOD_PATTERN ...	/* pattern only, no numerical values */
#define CHOLMOD_REAL ...		/* a real matrix */
#define CHOLMOD_COMPLEX ...	/* a complex matrix (ANSI C99 compatible) */
#define CHOLMOD_ZOMPLEX ...	/* a complex matrix (MATLAB compatible) */

/* -------------------------------------------------------------------------- */
/* cholmod_free_triplet:  free a triplet matrix */
/* -------------------------------------------------------------------------- */
int cholmod_l_free_triplet
(
    /* ---- in/out --- */
    cholmod_triplet **T,    /* triplet matrix to deallocate, NULL on output */
    /* --------------- */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_check_triplet:  check a sparse matrix in triplet form */
/* -------------------------------------------------------------------------- */
// This one is from "cholmod_check.h".
// Returns TRUE (1) if successful, or FALSE (0) otherwise.
int cholmod_l_check_triplet
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* triplet matrix to check */
    /* --------------- */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_print_triplet:  print a triplet matrix */
/* -------------------------------------------------------------------------- */
int cholmod_l_print_triplet
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* triplet matrix to print */
    const char *name,	/* printed name of triplet matrix */
    /* --------------- */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_sparse_to_triplet:  create a triplet matrix copy of a sparse matrix*/
/* -------------------------------------------------------------------------- */
cholmod_triplet *cholmod_l_sparse_to_triplet
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to copy */
    /* --------------- */
    cholmod_common *Common
) ;

/* -------------------------------------------------------------------------- */
/* cholmod_triplet_to_sparse:  create a sparse matrix copy of a triplet matrix*/
/* -------------------------------------------------------------------------- */
cholmod_sparse *cholmod_l_triplet_to_sparse
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* matrix to copy */
    size_t nzmax,	/* allocate at least this much space in output matrix */
    /* --------------- */
    cholmod_common *Common
) ;

// We need to call: cholmod_l_free( A.ncols, sizeof (SuiteSparse_long), E, cc ) ;
// UPDATE: Apparently we don't.
void *cholmod_l_free	/* always returns NULL */
(
    /* ---- input ---- */
    size_t n,		/* number of items */
    size_t size,	/* size of each item */
    /* ---- in/out --- */
    void *p,		/* block of memory to free */
    /* --------------- */
    cholmod_common *Common
) ;

/* cs_spsolve
int cs_spsolve
{
*/



/* [Z,R,E] = rz(A), returning Z, R, and E */
SuiteSparse_long SuiteSparseQR_C /* returns ???, (-1) if failure */
(
    /* inputs: */
    int ordering,               /* all, except 3:given treated as 0:fixed */
    double tol,                 /* columns with 2-norm <= tol treated as 0 */
    SuiteSparse_long econ,      /* e = max(min(m,econ),rank(A)) */
    int getCTX,                 /* 0:Z=C, 1:Z=c', 2: z=X */
    cholmod_sparse *A,          /* m-by-n sparse matrix to factorize */
    cholmod_sparse *Bsparse,
    cholmod_dense *Bdense,
    /* outputs: */
    cholmod_sparse **Zsparse,   /* m-by-e sparse matrix */
    cholmod_dense **Zdense,
    cholmod_sparse **R,         /* e-by-n sparse matrix */
    SuiteSparse_long **E,       /* size n column perm, NULL if identity */
    cholmod_sparse **H,         /* m-by-nh Householder vectors */
    SuiteSparse_long **HPinv,   /* size m row permutation */
    cholmod_dense **HTau,       /* 1-by-nh Householder coefficients */
    cholmod_common *cc          /* workspace and parameters */
) ;


/* [Q,R,E] = qr(A), returning Q as a sparse matrix */
SuiteSparse_long SuiteSparseQR_C_QR /* returns rank(A) est., (-1) if failure */
(
    /* inputs: */
    int ordering,               /* all, except 3:given treated as 0:fixed */
    double tol,                 /* columns with 2-norm <= tol treated as 0 */
    SuiteSparse_long econ,      /* e = max(min(m,econ),rank(A)) */
    cholmod_sparse *A,          /* m-by-n sparse matrix to factorize */
    /* outputs: */
    cholmod_sparse **Q,         /* m-by-e sparse matrix */
    cholmod_sparse **R,         /* e-by-n sparse matrix */
    SuiteSparse_long **E,       /* size n column perm, NULL if identity */
    cholmod_common *cc          /* workspace and parameters */
) ;


/* X = A\B where B is dense */
cholmod_dense *SuiteSparseQR_C_backslash    /* returns X, NULL if failure */
(
    int ordering,               /* all, except 3:given treated as 0:fixed */
    double tol,                 /* columns with 2-norm <= tol treated as 0 */
    cholmod_sparse *A,          /* m-by-n sparse matrix */
    cholmod_dense  *B,          /* m-by-k */
    cholmod_common *cc          /* workspace and parameters */
) ;

 /* X = A\B where B is sparse */
cholmod_sparse *SuiteSparseQR_C_backslash_sparse   /* returns X, or NULL */
(
    /* inputs: */
    int ordering,               /* all, except 3:given treated as 0:fixed */
    double tol,                 /* columns with 2-norm <= tol treated as 0 */
    cholmod_sparse *A,          /* m-by-n sparse matrix */
    cholmod_sparse *B,          /* m-by-k */
    cholmod_common *cc          /* workspace and parameters */
) ;

#define SPQR_ORDERING_FIXED ...
#define SPQR_ORDERING_NATURAL ...
#define SPQR_ORDERING_COLAMD ...
#define SPQR_ORDERING_GIVEN ...       /* only used for C/C++ interface */
#define SPQR_ORDERING_CHOLMOD ...     /* CHOLMOD best-effort (COLAMD, METIS,...)*/
#define SPQR_ORDERING_AMD ...         /* AMD(A'*A) */
#define SPQR_ORDERING_METIS ...       /* metis(A'*A) */
#define SPQR_ORDERING_DEFAULT ...     /* SuiteSparseQR default ordering */
#define SPQR_ORDERING_BEST ...        /* try COLAMD, AMD, and METIS; pick best */
#define SPQR_ORDERING_BESTAMD ...     /* try COLAMD and AMD; pick best */

#define SPQR_DEFAULT_TOL ...       /* if tol <= -2, the default tol is used */
#define SPQR_NO_TOL ...            /* if -2 < tol < 0, then no tol is used */
""")

def main():
    ## Two dirnames because ffibuilder.set_source()
    ## passes the module name: "sparseqr._sparseqr"
    ffibuilder.compile( verbose = True, tmpdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) )

if __name__ == "__main__":
    main()
