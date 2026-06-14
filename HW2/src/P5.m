load("cgn_illustration.mat");

[~,~,X_gmres,~,timer_gmres] = gmres_arnoldi(A, b, 70);
[X_cgn, timer_cgn] = cgn(A, b, 70);

figure
error_gmres = vecnorm(A*X_gmres-b);
error_cgn = vecnorm(A*X_cgn-b);

semilogy(error_gmres, "k--")
hold on
semilogy(error_cgn, "r-")
hold off
legend("GMRES", "CGN")
xlabel("Iteration")
ylabel("norm(Ax-b)")
xlim([1, 50])

figure
%[~, indices_gmres] = arrayfun(@(x) min(abs(timer_gmres - x)), target_values);
%[~, indices_cgn] = arrayfun(@(x) min(abs(timer_cgn - x)), target_values);

%ans_gmres = X_gmres(:,indices_gmres);
%ans_cgn = X_cgn(:,indices_cgn);

semilogy(timer_gmres, vecnorm(A*X_gmres-b), "k*", "DisplayName", "GMRES")
hold on
semilogy(timer_gmres, vecnorm(A*X_gmres-b), "k-")
semilogy(timer_cgn, vecnorm(A*X_cgn-b), "r*", "DisplayName", "CGN")
semilogy(timer_cgn, vecnorm(A*X_cgn-b), "r-")
hold off
legend()
xlabel("CPU-time (s)")
ylabel("norm(Ax-b)")
xlim([0, 6])



