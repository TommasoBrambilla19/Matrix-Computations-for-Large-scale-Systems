# Homework 2: Iterative Solvers for Linear Systems

This homework shifts focus to solving massive, sparse linear systems ($Ax=b$) using iterative Krylov subspace solvers where direct inversion is impossible.

### Implementations
All scripts are located in the `src/` directory:
* **`gmres_arnoldi.m`**: A custom Generalized Minimal Residual Method (GMRES) solver leveraging the Arnoldi iterations built in Homework 1.
* **`cg.m`**: Implementation of the Conjugate Gradient (CG) algorithm, optimized for symmetric positive-definite (SPD) matrices.
* **`cgn.m`**: Extension of the CG approach applied implicitly to the normal equations ($A^T A x = A^T b$) to handle non-symmetric matrices.
* **`cg_lsq.m`**: Adaptation of the Conjugate Gradient method to iteratively converge on optimal least squares solutions for overdetermined systems.
* **`P1.m`, `P2.m`, `P5.m`**: Analysis and execution scripts comparing the convergence rates, residual norms, and CPU times of the solvers above.