# Discrete Laplacian Distribution

This repository contains a header only library to sample from the discrete laplacian distribution in C++. The distribution is defined in "A discrete analogue of the Laplace distribution" by Seidu Inusaha and Tomasz J. Kozubowski. 

The interface follows the C++ STL interface for distributions. There also exists a test suite for confirming that the distribution of sampled values behaves as the discrete Laplacian distribution.

## Requirements

The library itself does not have any dependencies outside of the STL, but it requires C++20. If you need to use a lower standard, please open an issue and I will provide a fix.

The tests depend on Catch2.
