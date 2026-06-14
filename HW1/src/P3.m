%% part B
clc
rand('seed', 0);
nn = 10;
A = gallery('wathen', nn, nn);
b = ones(length(A(:,1)),1);

ylim([0 500]);
hold on
for m=1:40
    disp(m)
    K = Krylov(A,b,m);

    B = (K'*K)\(K'*A*K);

    ans2 = eig(B);
    scatter(m*ones(m,1), ans2, 'rx');


    [Q, H] = arnoldi(A,b,m);
    ansArn = eig(H(1:m,1:m));
    scatter(m*ones(m,1), ansArn, 'k.');
end

function out=Krylov(A,b,m)
    out = zeros(length(b), m);
    out(:,1) = b;
    for i = 2:m
        out(:,i) = A^(i)*b;
    end
end