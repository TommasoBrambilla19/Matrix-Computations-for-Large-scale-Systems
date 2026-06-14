# Homework 1: Eigenvalue Algorithms & Krylov Subspaces

This section explores the fundamental algorithms for eigenvalue extraction and builds the foundation for Krylov subspace methods. 

### Implementations
All scripts are located in the `src/` directory:
* **`P1.m`**: Implementation of the Power Method and Rayleigh Quotient Iteration, analyzing theoretical vs. practical convergence rates.
* **`P2.m`**: Benchmarking orthogonalization stability. Compares Classical, Modified, Double, and Triple Gram-Schmidt procedures by measuring orthogonality loss due to floating-point round-off errors.
* **`P3.m`**: Demonstrates the theoretical equivalence and practical numerical differences between primitive Krylov orthogonalization and Arnoldi factorization.
* **`P4.m`**: Implements the Shift-and-Invert Arnoldi Method to accelerate convergence for clustered and interior eigenvalues.
* **`arnoldi.m`**: A custom, robust implementation of the Arnoldi iteration used as a foundation for subsequent solvers.