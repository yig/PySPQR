# -*- coding: utf-8 -*-
#
"""setuptools-based setup.py for PySPQR.

Supports Python 2.7 and 3.4.

Usage as usual with setuptools:
    python setup.py build
    python setup.py sdist
    python setup.py bdist_wheel --universal  # use --universal for projects that work on both py2/py3 as-is
    python setup.py install

For details, see
    http://setuptools.readthedocs.io/en/latest/setuptools.html#command-reference
or
    python setup.py --help
    python setup.py --help-commands
    python setup.py --help bdist_wheel  # or any command
"""

from __future__ import division, print_function, absolute_import

try:
    # Python 3
    MyFileNotFoundError = FileNotFoundError
except:  # FileNotFoundError does not exist in Python 2.7
    # Python 2.7
    # - open() raises IOError
    # - remove() (not currently used here) raises OSError
    MyFileNotFoundError = (IOError, OSError)

#########################################################
# General config
#########################################################

# Name of the top-level package of your library.
#
# This is also the top level of its source tree, relative to the top-level project directory setup.py resides in.
#
libname="sparseqr"

# Short description for package list on PyPI
#
SHORTDESC="Python wrapper for SuiteSparseQR"

# Long description for package homepage on PyPI
#
DESC="""This module wraps the SuiteSparse QR decomposition and QR-based sparse linear solver functions for use with SciPy.

This is Matlab's sparse `[Q,R,E] = qr()`.

The solvers work also for overdetermined linear systems, making them useful for solving linear least-squares problems.

Supports Python 2.7 and 3.4.
"""

# Set up data files for packaging.
#
# Directories (relative to the top-level directory where setup.py resides) in which to look for data files.
#
# data_files, when used for this purpose, doesn't play nice with Mac OS (setuptools prefix set to /usr/local); we use MANIFEST.in instead.
#
# The MANIFEST.in solution makes README.md and test/test.py to be included in the sdist, but not the bdist (as they are not inside any package).
# From the user's perspective, this is probably fine.
#
# For some documentation, see
#   http://blog.cykerway.com/posts/2016/10/14/install-package-data-with-setuptools.html
#   https://stackoverflow.com/questions/24727709/i-dont-understand-python-manifest-in
#
# Note especially the comment
#   https://stackoverflow.com/questions/24727709/i-dont-understand-python-manifest-in#comment46482024_24727824
# which states that
#    To head off the inevitable package_data and data_files recommendations, which are out of scope, I'll continue.
#    package_data lists file that get installed with your package into dist-packages/yourpackage which would have been
#    skipped because the don't have a *.py name. data_files lists files that get installed outside of your package.
#    Each entry specifies a target path that is prefixed with sys.prefix if it is relative or created directly
#    (permissions permitting) if it begins with a /. -Bruno Bronosky Mar 18 '15 at 16:55 
#
datadirs  = []

# File extensions to be considered as data files. (Literal, no wildcards.)
dataexts  = (".py", ".ipynb",  ".sh",  ".lyx", ".tex", ".txt", ".pdf")

# Standard documentation to detect (and package if it exists).
#
#standard_docs     = ["README", "LICENSE", "TODO", "CHANGELOG", "AUTHORS"]  # just the basename without file extension
standard_docs     = []  # same here, don't use data_files
standard_doc_exts = [".md", ".rst", ".txt", ""]  # commonly .md for GitHub projects, but other projects may use .rst or .txt (or even blank).


#########################################################
# Init
#########################################################

# check for Python 2.7 or later
# http://stackoverflow.com/questions/19534896/enforcing-python-version-in-setup-py
import sys
if sys.version_info < (2,7):
    sys.exit('Sorry, Python < 2.7 is not supported')

import os

from setuptools import setup


# Gather user-defined data files
#
# http://stackoverflow.com/questions/13628979/setuptools-how-to-make-package-contain-extra-data-folder-and-all-folders-inside
#
datafiles = []
getext = lambda filename: os.path.splitext(filename)[1]
for datadir in datadirs:
    datafiles.extend( [(root, [os.path.join(root, f) for f in files if getext(f) in dataexts])
                       for root, dirs, files in os.walk(datadir)] )


# Add standard documentation (README et al.), if any, to data files
#
detected_docs = []
for docname in standard_docs:
    for ext in standard_doc_exts:
        filename = "".join( (docname, ext) )  # relative to the directory in which setup.py resides
        if os.path.isfile(filename):
            detected_docs.append(filename)
datafiles.append( ('.', detected_docs) )


# Extract __version__ from the package __init__.py
# (since it's not a good idea to actually run __init__.py during the build process).
#
# http://stackoverflow.com/questions/2058802/how-can-i-get-the-version-defined-in-setup-py-setuptools-in-my-package
#
import ast
init_py_path = os.path.join(libname, '__init__.py')
version = '0.0.unknown'
try:
    with open(init_py_path) as f:
        for line in f:
            if line.startswith('__version__'):
                version = ast.parse(line).body[0].value.s
                break
        else:
            print( "WARNING: Version information not found in '%s', using placeholder '%s'" % (init_py_path, version), file=sys.stderr )
except MyFileNotFoundError:
    print( "WARNING: Could not find file '%s', using placeholder version information '%s'" % (init_py_path, version), file=sys.stderr )


#########################################################
# Call setup()
#########################################################

setup(
    name = "sparseqr",
    version = version,
    author = "Yotam Gingold and Juha Jeronen",
#    author_email = "juha.jeronen@tut.fi",
    url = "https://github.com/yig/PySPQR",

    description = SHORTDESC,
    long_description = DESC,

    license = "Public Domain CC0",

    # free-form text field; http://stackoverflow.com/questions/34994130/what-platforms-argument-to-setup-in-setup-py-does
    platforms = ["Linux"],

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
    cffi_modules = ["sparseqr/spqr_gen.py:ffibuilder"],
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

    # This flag is usually used with MANIFEST.in, but we don't need it (our README and test.py do not reside in any package and thus would not be installed in a bdist anyway; they are only to be included in the sdist)
#    include_package_data = True

    # Custom data files not inside a Python package
#    data_files = datafiles
)

