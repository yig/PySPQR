# -*- coding: utf-8 -*-
#
from __future__ import division, print_function, absolute_import

from setuptools import setup

setup(
    name = "sparseqr",
    version = "1.0.0",
    author = "Yotam Gingold and Juha Jeronen",
    url = "https://github.com/yig/PySPQR",

    description = "Python wrapper for SuiteSparseQR",
    long_description = """This module wraps the SuiteSparse QR decomposition and QR-based sparse linear solver functions for use with SciPy.

This is Matlab's sparse `[Q,R,E] = qr()`.

The solver works also for overdetermined linear systems, making it useful for solving linear least-squares problems.

Supports Python 2.7 and 3.4.
""",

    license = "Public Domain CC0",

    # free-form text field; http://stackoverflow.com/questions/34994130/what-platforms-argument-to-setup-in-setup-py-does
    platforms = ["any"],

    # See
    #    https://pypi.python.org/pypi?%3Aaction=list_classifiers
    #
    # for the standard classifiers.
    #
    classifiers = [ "Development Status :: 4 - Beta",
                    "Environment :: Console",
                    "Intended Audience :: Developers",
                    "Intended Audience :: Science/Research",
                    "License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication",
                    "Operating System :: POSIX :: Linux",
                    "Programming Language :: Python",
                    "Programming Language :: Python :: 2",
                    "Programming Language :: Python :: 2.7",
                    "Programming Language :: Python :: 3",
                    "Programming Language :: Python :: 3.4",
                    "Topic :: Scientific/Engineering",
                    "Topic :: Scientific/Engineering :: Mathematics",
                    "Topic :: Software Development :: Libraries",
                    "Topic :: Software Development :: Libraries :: Python Modules"
                  ],

    # See
    #    http://setuptools.readthedocs.io/en/latest/setuptools.html
    #
    setup_requires = ["cffi>=1.0.0"],
    cffi_modules = ["sparseqr/sparseqr_gen.py:ffibuilder"],
    install_requires = ["numpy", "scipy", "cffi>=1.0.0"],
    provides = ["sparseqr"],

    # keywords for PyPI (in case you upload your project)
    #
    # e.g. the keywords your project uses as topics on GitHub, minus "python" (if there)
    #
    keywords = ["suitesparse bindings wrapper scipy numpy qr-decomposition qr-factorisation sparse-matrix sparse-linear-system sparse-linear-solver"],

    # Declare packages so that  python -m setup build  will copy .py files (especially __init__.py).
    #
    # This **does not** automatically recurse into subpackages, so they must also be declared.
    #
    packages = ["sparseqr"],

    zip_safe = False  # includes a binary extension, not zip safe
)
