# Python wrapper for SuiteSparseQR

This module wraps the [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html)
decomposition function for use with [SciPy](http://www.scipy.org).
This is Matlab's sparse `[Q,R,E] = qr()`.
For some reason, no one ever wrapped that function of SuiteSparseQR for Python.

Also wrapped are the SuiteSparseQR solvers for ``A x = b`` for the cases with sparse `A` and dense or sparse `b`.
This is especially useful for solving sparse overdetermined linear systems in the least-squares sense.
Here `A` is of size m-by-n and `b` is m-by-k (storing `k` different right-hand side vectors, each considered separately).

# Usage

```python
import numpy
import scipy.sparse.linalg
import sparseqr

# QR decompose a sparse matrix M such that  Q R = M E
#
M = scipy.sparse.rand( 10, 10, density = 0.1 )
Q, R, E, rank = sparseqr.qr( M )
print( abs( Q*R - M*sparseqr.permutation_vector_to_matrix(E) ).sum() )  # should be approximately zero

# Solve many linear systems "M x = b for b in columns(B)"
#
B = scipy.sparse.rand( 10, 5, density = 0.1 )  # many RHS, sparse (could also have just one RHS with shape (10,))
x = sparseqr.solve( M, B, tolerance = 0 )

# Solve an overdetermined linear system  A x = b  in the least-squares sense
#
# The same routine also works for the usual non-overdetermined case.
#
A = scipy.sparse.rand( 20, 10, density = 0.1 )  # 20 equations, 10 unknowns
b = numpy.random.random(20)  # one RHS, dense, but could also have many (in shape (20,k))
x = sparseqr.solve( A, b, tolerance = 0 )

# Solve a linear system  M x = B  via QR decomposition
#
# This approach is slow due to the explicit construction of Q, but may be
# useful if a large number of systems need to be solved with the same M.
#
M = scipy.sparse.rand( 10, 10, density = 0.1 )
Q, R, E, rank = sparseqr.qr( M )
r = rank  # r could be min(M.shape) if M is full-rank

# The system is only solvable if the lower part of Q.T @ B is all zero:
print( "System is solvable if this is zero:", abs( (( Q.tocsc()[:,r:] ).T ).dot( B ) ).sum() )

# Systems with large non-square matrices can benefit from "economy" decomposition.
M = scipy.sparse.rand( 20, 5, density=0.1 )
B = scipy.sparse.rand( 20, 5, density = 0.1 )
Q, R, E, rank = sparseqr.qr( M )
print("Q shape:", Q.shape)  # Q shape: (20, 20)
print("R shape:", R.shape)  # R shape: (20, 5)
Q, R, E, rank = sparseqr.qr( M, economy=True )
print("Q shape:", Q.shape)  # Q shape: (20, 5)
print("R shape:", R.shape)  # R shape: (5, 5)

# Use CSC format for fast indexing of columns.
R = R.tocsc()[:r,:r]
Q = Q.tocsc()[:,:r]
QB = (Q.T).dot(B).tocsc()  # for best performance, spsolve() wants the RHS to be in CSC format.
result = scipy.sparse.linalg.spsolve(R, QB)

# Recover a solution (as a dense array):
x = numpy.zeros( ( M.shape[1], B.shape[1] ), dtype = result.dtype )
x[:r] = result.todense()
x[E] = x.copy()

# Recover a solution (as a sparse matrix):
x = scipy.sparse.vstack( ( result.tocoo(), scipy.sparse.coo_matrix( ( M.shape[1] - rank, B.shape[1] ), dtype = result.dtype ) ) )
x.row = E[ x.row ]
```

# Installation

## Via `pip`

From PyPI:

```bash
pip install sparseqr
```

From GitHub:

```bash
pip install git+https://github.com/yig/PySPQR.git
```

## Manually from GitHub

As user:

```bash
git clone https://github.com/yig/PySPQR.git
cd PySPQR
python setup.py install --user
```

As admin, change the last command to

```bash
sudo python setup.py install
```

## Directly

Copy the three `sparseqr/*.py` files next to your source code,
or leave them in their directory and call it as a module.


# Deploy

Change the version in:

    sparseqr/__init__.py
    setup.py 
    pyproject.toml 

Run:

    poetry build -f sdist
    poetry publish

# Tested on

 - Python 2.7, 3.4, 3.5, and 3.9.
 - Conda and not conda.
 - Mac OS X, Ubuntu Linux and Linux Mint.


# Dependencies

* [SciPy/NumPy](http://www.scipy.org)
* [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html) (macOS: `brew install suitesparse`; debian/ubuntu linux: `apt-get install libsuitesparse-dev`)
* [cffi](http://cffi.readthedocs.io/) (`pip install cffi`)

# License

Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
