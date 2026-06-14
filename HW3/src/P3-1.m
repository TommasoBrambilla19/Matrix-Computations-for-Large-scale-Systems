%% A
A = [1 4 4; 3 -1 3; -1 4 4];
f = @(x) sin(x);

schur_parlett(A, f);

sin(A);
%% B
A=rand(100,100)+1i*rand(100,100);
A=A/norm(A);

rep = 10;
N_vec = 1:500;
timers = zeros(length(N_vec), 2);
for r=1:rep
    for n=N_vec
        % Naive
        tic
        B=A; for i=1:n-1; B=B*A; end
        timers(n,1) = timers(n,1) + toc;

        % Schur-Parlett
        tic
        f = @(x) x^n;
        schur_parlett(A, f);
        timers(n,2) = timers(n,2) + toc;
    end
end

timers = timers/rep;
plot(timers(:,1), 'kx', 'DisplayName', 'Naive');
hold on;
plot(timers(:,2), 'ro', 'DisplayName', 'Schur-Parlett');
legend();