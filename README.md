# Discrete Laplacian Distribution

This repository contains a header only C++ library to sample from distributions used for differential privacy applications. Currently it includes:
- the discrete Laplacian distribution (defined in "A discrete analogue of the Laplace distribution" by Seidu Inusaha and Tomasz J. Kozubowski, using the algorithm outlined there)
- the discrete Gaussian distribution (using the algorithm outlined in "The Discrete Gaussian for Differential Privacy" by Canonne et al.

The interface follows the C++ STL interface for distributions, and implements the RandomNumberGenerator trait. There also exists a test suite for confirming that the distribution of sampled values behaves as they should, and some benchmarking code.

## Requirements

The library itself does not have any dependencies outside of the STL, but it requires C++20. If you need to use a lower standard, please open an issue and I will provide a fix.

The testing and benchmarking depends on [Catch2](https://github.com/catchorg/Catch2).
