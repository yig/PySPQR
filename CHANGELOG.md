# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [v1.5] - 2026-02-03
### Added
- Added complex matrix support in `cholmodsparse2scipy` and `scipy2cholmodsparse`.
- Modern test cases that can be run with `uv run --extra test pytest` and run automatically with continuous integration on GitHub.

## [v1.4.1] - 2025-02-10
### Fixed
- An import statement got lost in the 1.4 update.

## [v1.4] - 2025-01-28
### Fixed
- Modernized the build system (`pyproject.toml`).
- Changed the way suite-sparse is found to be more robust.

## [v1.3] - 2025-01-09
### Added
- Bindings for `qr_factorize` and `qmult` (thanks to jkrokowski)

### Fixed
- Compatibility with more environments (more search paths, newer numpy, setuptools dependency)
- Readme example uses `spsolve_triangular`.

## [v1.2.1] - 2023-04-12
### Fixed
- Fixed a memory leak in `qr()` and `rz()`.
### Changed
- Bumped minimal Python version to 3.8.
- `rz()` is called by the test script. Its output is ignored.

## [v1.2] - 2022-05-27
### Added
- Added support for partial "economy" decompositions. (Christoph Hansknecht <c.hansknecht@tu-braunschweig.de>): 'The "economy" option can be used in SPQR to compute a QR factorization of a (m x n) matrix with m < n consisting of blocks Q_1, and Q_2, where Q_1 has as shape of (m x n) and Q_2 of (m x k - n). For k = n we get the reduced form, for k = m the full one. For k in between m and n, SPQR yields a block that spans part of the kernel of A. This patch adds this functionality to PySPQR.'
- Added support for macOS on arm64.

## [v1.1.2] - 2021-08-09
### Added
- Added rz recomposition (thanks to Ben Smith <bsmith@apl.washington.edu>)
- Added support for "economy" decomposition. (Jeffrey Bouas <ignirtoq@gmail.com>)
### Changed
- Supports conda environments (thanks to Ben Smith <bsmith@apl.washington.edu> and Sterling Baird <sterling.baird@icloud.com>)

## [v1.0.0] - 2017-08-31
### Added
 - Installation and packaging using `setuptools`
### Changed
 - Rename module `spqr` to `sparseqr`
 - Clean up public API: `qr`, `solve`, `permutation_vector_to_matrix`

## [v1.0.0] - 2017-08-31
### Added
 - Installation and packaging using `setuptools` (thanks to Juha Jeronen <juha.jeronen@tut.fi>)
### Changed
 - Rename module `spqr` to `sparseqr`
 - Clean up public API: `qr`, `solve`, `permutation_vector_to_matrix`
