n=100;
e = ones(n,1);
A = -spdiags([e -6*e e], -1:1, n, n);
b = kron(ones(n/2,1),[0;1])+ones(n,1);

X_cg = cg(A, b, 30);
X_lsqcg = cg_lsq(A, b, 30);
[~,~,X_gmres,~,~] = gmres_arnoldi(A, b, 30);

figure
semilogy(vecnorm(A*X_cg-b), "b*")
hold on
semilogy(vecnorm(A*X_gmres-b), "ro")
semilogy(vecnorm(A*X_lsqcg-b), "k+")
hold off

legend("Method 1: CG", "Method 2: GMRES", "Method 3: CG with LSQ")
xlabel("Iteration")
ylabel("norm(Ax-b)")