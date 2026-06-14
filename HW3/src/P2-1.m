% A B C
rep = 10;
Ms = [10,100,200,300,400];
times = zeros(length(Ms),2);
for r = 1:rep
    for i = 1:length(Ms)
        A = alpha_example(1,Ms(i));
        tic;
        H1 = naive_hessenberg_red(A);
        times(i,2) = times(i,2) + toc;
        tic;
        H2 = hessenberg_red(A);
        times(i,1) = times(i,1) + toc;
    end
end

times = times/rep;

% D
format long
sigms = [0, 1];
epsilons = [0.4, 10.^(-(1:10)), 0];
values = zeros(length(epsilons), length(sigms));

for s = 1:length(sigms)
    for i = 1:length(epsilons)
        A = [3, 2; epsilons(i), 1];
        [eigvals, last_A, k] = QR_shift(A, sigms(s), 1, 1e-10);
        values(i,s) = abs(last_A(2,1));
    end
end

