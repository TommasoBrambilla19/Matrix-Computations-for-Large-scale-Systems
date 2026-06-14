format long
epsilons = logspace(-10, -1, 1000);
error_vec = zeros(1, length(epsilons));
rep = 100;

for r=1:rep
    for i=1:length(epsilons)
        A = [pi 1; 0 epsilons(i)+pi];

        % Linearized solution
        l1 = A(1,1); % matrix is upper triangular
        l2 = A(2,2); 
        g1 = exp(l1); % function of interest
        g2 = exp(l2);

        alpha = (g1*l2 - g2*l1)/(l2-l1);
        beta = (g2 - g1)/(l2-l1);
        p = @(z) alpha*eye(size(z)) + beta*z;

        fA = p(A);
        
        % Analytical solution
        [V,D]=eig(A);
        F=diag(exp(diag(D)));
        F=V*F*inv(V); 

        %Same as expm(A)
        %F = expm(A);

        error_vec(i) = error_vec(i) + norm(F-fA);
    end
end
error_vec = error_vec/rep;

loglog(epsilons, error_vec)
xlabel('$\epsilon$', 'Interpreter', 'latex');
ylabel('$||f(A)-F||$', 'Interpreter', 'latex');
xlim([1e-10, 1e-1]);