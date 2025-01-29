# -*- coding: utf-8 -*-
#
from __future__ import division, print_function, absolute_import

from setuptools import setup #, dist
import os

setup(
    # See
    #    http://setuptools.readthedocs.io/en/latest/setuptools.html
    #
    setup_requires = ["cffi>=1.0.0"],
    cffi_modules = ["sparseqr/sparseqr_gen.py:ffibuilder"],
    install_requires = ["cffi>=1.0.0"],
)
