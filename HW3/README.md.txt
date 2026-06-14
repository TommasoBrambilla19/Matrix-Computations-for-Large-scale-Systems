# Homework 3: The QR Algorithm & Matrix Functions

The final section explores advanced eigenvalue solvers and the highly complex task of evaluating functions of matrices.

### Implementations
All scripts are located in the `src/` directory:
* **`QR.m` & `QR_shift.m`**: Implementations of the standard QR algorithm and the shifted QR algorithm (using Givens rotations) to significantly accelerate eigenvalue convergence.
* **`naive_hessenberg_red.m` & `hessenberg_red.m`**: Algorithms for preprocessing matrices into upper Hessenberg form to optimize the QR algorithm, demonstrating the difference between naive and optimized memory access.
* **`schur_parlett.m`**: A solver for evaluating general scalar functions on matrices via the complex Schur decomposition.
* **`alpha_example.m`**: A matrix generator used to test algorithm bounds and behavior.
* **`P1.m` - `P5.m`**: Analytical scripts. Highlights include evaluating the matrix exponential ($e^A$) via Taylor series and the highly efficient "scaling and squaring" technique.