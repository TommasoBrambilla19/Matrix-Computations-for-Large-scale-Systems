function [x_lsq_stars] = cg_lsq(A, b, m)
    [~, num_cols] = size(A);
    x_m0 = zeros(num_cols, 1);
    r_m0 = b;
    p_m0 = r_m0;
    x_stars = zeros(num_cols, m);
    x_lsq_stars = zeros(num_cols, m);

    for i = 1:m
        a_m1 = r_m0'*r_m0/(p_m0'*A*p_m0);
        x_m1 = x_m0 + a_m1*p_m0;
        r_m1 = r_m0 - a_m1*A*p_m0;
        b_m1 = r_m1'*r_m1/(r_m0'*r_m0);
        p_m1 = r_m1 + b_m1*p_m0;

        p_m0 = p_m1;
        r_m0 = r_m1;
        x_m0 = x_m1;

        x_stars(:,i) = x_m0;
        X = x_stars(:,1:i);
        AX = A*X;
        x_lsq_stars(:,i) = X*((AX'*AX)\(AX'*b));
    end
end