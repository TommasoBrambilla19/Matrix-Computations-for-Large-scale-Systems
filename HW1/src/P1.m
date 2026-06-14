clear all
format long
A = [1 2 3; 
     2 0 2; 
     3 2 9];

x0 = [1; 0; -1];

[V, D] = eig(A);

%% Part 1 A
[lambda, w, l_vecs] = power_method(A, x0, 10);
[sorted, index] = sort(abs(diag(D)), 'descend');
ref = V(:,index(1));

err_vec = min([vecnorm(l_vecs - ref); vecnorm(l_vecs + ref)]);

semilogy(1:length(err_vec), err_vec);
hold on
semilogy(1:length(err_vec), ((sorted(2))/sorted(1)).^(1:length(err_vec)), "k--")
hold off

%% Part 1 B

[lambda_b, w_b, l_vals_b, l_vecs_b] = rayleigh_quotient_iteration(A, x0, 5);
err_vals_b = abs(l_vals_b(1:end-1)-l_vals_b(end));

err_vecs_b = zeros(1, length(l_vecs_b));
for i = 1:(length(l_vecs_b)-1)
    err_vecs_b(i)=min(norm(l_vecs_b(:,i)-l_vecs_b(:,end)), norm(l_vecs_b(:,i)+l_vecs_b(:,end)));
end
figure
disp(log(err_vecs_b(3:end)./err_vecs_b(2:end-1))./log(err_vecs_b(2:end-1)./err_vecs_b(1:end-2)))
semilogy(1:length(err_vals_b), err_vals_b, DisplayName='Original lambdas');
hold on

semilogy(1:length(err_vecs_b), err_vecs_b, 'DisplayName', 'Reference ^1');
semilogy(1:length(err_vecs_b), err_vecs_b.^2, 'DisplayName', 'Reference ^2');
semilogy(1:length(err_vecs_b), err_vecs_b.^3, 'DisplayName', 'Reference ^3');

legend(Location='southwest');
ylabel('Error')
xlabel('Iteration')




%% Part 1 C
B = [1 2 4;
     2 0 2; 
     3 2 9];

[lambda_c, w_c, l_vals_c, l_vecs_c] = rayleigh_quotient_iteration(B, x0, 6);
err_vals_c = abs(l_vals_c(1:end-1)-l_vals_c(end));
err_vecs_c = zeros(1, length(l_vecs_c));
for i = 1:(length(l_vecs_c)-1)
    err_vecs_c(i)=min(norm(l_vecs_c(:,i)-l_vecs_c(:,end)), norm(l_vecs_c(:,i)+l_vecs_c(:,end)));
end
figure
disp(log(err_vecs_c(3:end)./err_vecs_c(2:end-1))./log(err_vecs_c(2:end-1)./err_vecs_c(1:end-2)))
semilogy(1:length(err_vals_c), err_vals_c, DisplayName='Modified lambdas');
hold on

semilogy(1:length(err_vecs_c), err_vecs_c, 'DisplayName', 'Reference ^1');
semilogy(1:length(err_vecs_c), err_vecs_c.^2, 'DisplayName', 'Reference ^2');
semilogy(1:length(err_vecs_c), err_vecs_c.^3, 'DisplayName', 'Reference ^3');

legend(Location='southwest');
ylabel('Error')
xlabel('Iteration')

%disp(err_vals_c)
%disp(err_vals_b)







function [lambda, w, l_vecs] = power_method(A, v0, max_iter)
    v0 = v0/norm(v0);
    l_vecs = zeros(length(A), max_iter);
    for i = 1:max_iter
        w = A*v0;
        v0 = w/norm(w);
        lambda = v0'*A*v0;

        l_vecs(:,i) = v0;
    end
end

function [lambda, w, l_vals, l_vecs] = rayleigh_quotient_iteration(A, v0, max_iter)
    v0 = v0/norm(v0);
    lambda = v0'*A*v0; % Rayleigh quotient

    l_vals = zeros(1, max_iter);
    l_vecs = zeros(3, max_iter);
    l_vals(1) = lambda;
    l_vecs(:,1) = v0;

    for i = 2:max_iter
        w = (A - lambda*eye(size(A)))\v0;
        v0 = w/norm(w);
        lambda = v0'*A*v0;

        l_vals(i) = lambda;
        l_vecs(:,i) = v0;
    end
end