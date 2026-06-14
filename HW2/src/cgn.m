function [x_stars, timer] = cgn(A, b, m)
    [~, num_cols] = size(A);
    x_m0 = zeros(num_cols, 1);
    r_m0 = A'*b;
    p_m0 = r_m0;
    x_stars = zeros(num_cols, m);
    timer = zeros(m,1);
    
    tic
    x_stars(:,1) = x_m0;
    for i = 1:m
        ATP = A'*(A*p_m0);
        a_m1 = r_m0'*r_m0/(p_m0'*ATP);
        x_m1 = x_m0 + a_m1*p_m0;
        r_m1 = r_m0 - a_m1*ATP;
        b_m1 = r_m1'*r_m1/(r_m0'*r_m0);
        p_m1 = r_m1 + b_m1*p_m0;

        p_m0 = p_m1;
        r_m0 = r_m1;
        x_m0 = x_m1;

        x_stars(:,i) = x_m1;
        timer(i) = toc;
    end
end