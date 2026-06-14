function [eigvals, last_A, i] = QR(A, maxIter, tol)
    [Q, R] = qr(A);
    for i = 1:maxIter
        [Q, R] = qr(A);
        A = R*Q;
        if max(max(abs(tril(A,-1)))) < tol || i == maxIter
            eigvals = diag(A);
            last_A = A;
            break;
        end
    end
end