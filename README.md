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
print( "Should be approximately zero:", abs( Q*R - M*sparseqr.permutation_vector_to_matrix(E) ).sum() ) 

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
## Call `rz()`:
sparseqr.rz( A, b, tolerance = 0 )

# Solve a linear system  M x = B  via QR decomposition
#
# This approach is slow due to the explicit construction of Q, but may be
# useful if a large number of systems need to be solved with the same M.
#
M = scipy.sparse.rand( 10, 10, density = 0.1 )
Q, R, E, rank = sparseqr.qr( M )
r = rank  # r could be min(M.shape) if M is full-rank

# The system is only solvable if the lower part of Q.T @ B is all zero:
print( "System is solvable if this is zero (unlikely for a random matrix):", abs( (( Q.tocsc()[:,r:] ).T ).dot( B ) ).sum() )

# Systems with large non-square matrices can benefit from "economy" decomposition.
M = scipy.sparse.rand( 20, 5, density=0.1 )
B = scipy.sparse.rand( 20, 5, density = 0.1 )
Q, R, E, rank = sparseqr.qr( M )
print("Q shape (should be 20x20):", Q.shape)
print("R shape (should be 20x5):", R.shape)
Q, R, E, rank = sparseqr.qr( M, economy=True )
print("Q shape (should be 20x5):", Q.shape)
print("R shape (should be 5x5):", R.shape)


R = R.tocsr()[:r,:r] #for best performance, spsolve_triangular() wants the Matrix to be in CSR format.
Q = Q.tocsc()[:,:r] # Use CSC format for fast indexing of columns.
QB = (Q.T).dot(B).todense()  # spsolve_triangular() need the RHS in array format.
result = scipy.sparse.linalg.spsolve_triangular(R, QB, lower=False)

# Recover a solution (as a dense array):
x = numpy.zeros( ( M.shape[1], B.shape[1] ), dtype = result.dtype )
x[:r] = result
x[E] = x.copy()

# Recover a solution (as a sparse matrix):
x = scipy.sparse.vstack( ( result, scipy.sparse.coo_matrix( ( M.shape[1] - rank, B.shape[1] ), dtype = result.dtype ) ) )
x.row = E[ x.row ]
```

# Installation

Before installing this module, you must first install [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html). You can do that via conda (`conda install suitesparse`) or your system's package manager (macOS: `brew install suitesparse`; debian/ubuntu linux: `apt-get install libsuitesparse-dev`).

Now you are ready to install this module.

## Via `pip`

From PyPI:

```bash
pip install sparseqr
```

From GitHub:

```bash
pip install git+https://github.com/yig/PySPQR.git
```

## Directly

Copy the three `sparseqr/*.py` files next to your source code,
or leave them in their directory and call it as a module.


# Deploy

1. Change the version in `sparseqr/__init__.py`

2. Update `CHANGELOG.md`

3. Commit to git. Push to GitHub.

4. Run (in a clean repos, e.g., `git clone . clean; cd clean`):

    ```
    flit publish --format sdist
    ```

    Using [uv](https://docs.astral.sh/uv/) and [PyPI API tokens](https://pypi.org/help/#apitoken):

    ```
    FLIT_USERNAME=__token__ uv tool run --with flit flit publish --format sdist
    ```

We don't publish binary wheels, because it must be compiled against suite-sparse as a system dependency. We could publish a `none-any` wheel, which would cause compilation to happen the first time the module is imported rather than when it is installed. Is there a point to that?

# Known issues

`pip uninstall sparseqr` won't remove the generated libraries. It will list them with a warning.

# Tested

GitHub Continuous Integration (CI) tests:

 - Python 3.9, 3.10, 3.11, 3.12, 3.13, 3.14.
 - macOS, Ubuntu Linux, and Windows.
 - conda (Windows) and not conda (Linux/macOS). I tested conda on macOS manually at one point.

Test manually with:

```
python -m pytest
```

or

```
uv run --extra test pytest
```

# Dependencies

These are listed as dependencies and will be installed automatically:

* [SciPy/NumPy](http://www.scipy.org)
* [cffi](http://cffi.readthedocs.io/)

These must be installed manually:

* [SuiteSparseQR](http://faculty.cse.tamu.edu/davis/suitesparse.html) (macOS: `brew install suitesparse`; debian/ubuntu linux: `apt-get install libsuitesparse-dev`)

# License

Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
