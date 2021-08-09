# -*- coding: utf-8 -*-
#
"""Python bindings for SuiteSparse QR routines.

Exported functions:
    sparseqr.qr                              QR decompose sparse matrix
    sparseqr.solve                           solve linear system, LHS sparse
    sparseqr.permutation_vector_to_matrix    utility for conversion

The solver works also for overdetermined linear systems,
making it useful for solving linear least-squares problems.

In solve(), the RHS can be dense or sparse.

See the docstrings of the individual functions for details.
"""

from __future__ import absolute_import

__version__ = '1.1.1'

# import the important things into the package's top-level namespace.
from .sparseqr import qr, rz, solve, permutation_vector_to_matrix

