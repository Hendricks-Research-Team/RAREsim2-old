# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## \[Unreleased]

### Added

* Prob mode now prints inferred probability/MAC bins

### Changed

* Sim now prints a single allele frequency summary table for all pruning modes
* Updated standard output formatting for sim and calc functions

### Fixed

* Haplotype parsing so trailing whitespace no longer creates an extra all-zero output column
* Probabilistic pruning to prune whole variant rows instead of individual alleles
* Intronic rows removed from example input legend and haplotype files

## \[3.0.2] - 2026-02-18

### Added

* Z flag support and pruned-variants-file for -prob mode
* New example files with monomorphic variants removed

### Changed

* Reorganized the example directory
* Updated documentation based on new examples

### Removed

* rearch\_example.zip

## \[3.0.0] - 2025-10-18

### Added

* Initial project setup with modern Python packaging
* Sparse matrix implementation for efficient storage of genetic data
* Command-line interface for common operations
* Comprehensive test suite
* Documentation and contribution guidelines

### Changed

* Migrated from legacy setup.py to pyproject.toml
* Improved code organization and structure
* Enhanced error handling and validation

### Fixed

* Various bug fixes and performance improvements

### Removed

* C/Cython dependencies



