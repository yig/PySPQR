# -*- coding: utf-8 -*-
#
"""Python bindings for SuiteSparse QR routines.

Exported functions:
    sparseqr.qr                    QR decomposition
    sparseqr.solve                 solve linear system; LHS sparse, RHS dense or sparse
    sparseqr.permutation_from_E    convert permutation vector to permutation matrix

The solver works also for overdetermined linear systems,
making it useful for solving linear least-squares problems.

See the docstrings of the individual functions for details.
"""

from __future__ import absolute_import

__version__ = '1.0.0'

# import the important things into the package's top-level namespace.
from .sparseqr import qr, solve, permutation_from_E

