function [eigvals, last_A, k] = QR_shift(A, shift, maxIter, tol) % Assumed A is hessenberg already
    n = length(A);
    G(1:n-1)={zeros(2)};
    shifter = shift*eye(n);

    for k = 1:maxIter % QR
        H = A-shifter;
        for i = 1:n-1
            [G{i}, ~] = planerot(H(i:i+1,i));
            H(i:i+1,i:n) = G{i}*H(i:i+1,i:n);
        end

        for i = 1:n-1
            H(1:i+1,i:i+1) = H(1:i+1,i:i+1)*G{i}';
        end

        if max(max(abs(tril(A,-1)))) < tol || k == maxIter
            eigvals = diag(A);
            last_A = H+shifter;
            break;
        end
        A = H+shifter;
    end

end