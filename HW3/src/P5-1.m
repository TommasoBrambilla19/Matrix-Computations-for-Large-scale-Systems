n=256; randn('seed',0);
A=randn(n,n);
A=A/norm(A); 
A=A-3*eye(n);

ref = expm(A);
error_mat = zeros(7, 8);
time_mat = zeros(7, 8);
rep = 10;

for r=1:rep
    for m=1:7
        for j=0:7
            tic
            i = 2^j;
            P_jm = Taylor(A/(i), m)^(i);
            error_mat(m,j+1) = error_mat(m,j+1) + log10(norm(P_jm-ref));
            time_mat(m,j+1) = time_mat(m,j+1) + toc;
        end
    end
end

error_mat_j0 = zeros(40, 1);
time_mat_j0 = zeros(40, 1);
for r=1:rep
    for m=1:40
        j=0;
        tic
        i = 2^j;
        P_jm = Taylor(A/(i), m)^(i);
        error_mat_j0(m,j+1) = error_mat_j0(m,j+1) + log10(norm(P_jm-ref));
        time_mat_j0(m,j+1) = time_mat_j0(m,j+1) + toc;
    end
end
error_mat_j0 = error_mat_j0/rep;
time_mat_j0 = time_mat_j0/rep;

error_mat = error_mat/rep;
time_mat = time_mat/rep;

function out = Taylor(A,m)
    out=eye(size(A));
    B=1;
    for k=1:m
        B=B*A;
        out = out + B/(factorial(k));
    end
end