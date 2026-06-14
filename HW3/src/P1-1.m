n = 1000;
alpha = linspace(1e1, 1e5, n);
iter_vec = zeros(n, 1);
theory_vec = zeros(n, 1);

for i = 1:n
    A = alpha_example(alpha(i));
    [~, ~, iter] = QR(A, 1e5, 1e-10);
    iter_vec(i) = iter;

    % Theoretical reasoning
    eigvals = sort(eig(A), "descend");
    dom_term = errfun(eigvals./eigvals');
    theory_vec(i) = log(1e-10)/log(dom_term);
end



semilogx(alpha, iter_vec, 'r-','DisplayName', 'Number of iterations');
hold on;
semilogx(alpha, theory_vec, 'k--','DisplayName', 'Predicted number of iterations');

xlabel('$\alpha$', 'Interpreter', 'latex');
ylabel('Number of iterations', 'Interpreter', 'latex');
legend();
ylim([0, 3000]);

function out = errfun(A)
    out=max(max(abs(tril(A,-1))));
end