function [x_stars] = cg(A, b, m)
    [~, num_cols] = size(A);
    x_m0 = zeros(num_cols, 1);
    x_stars = zeros(num_cols, m);
    
    r_m0 = b;
    p_m0 = r_m0;
    
    for i = 1:m
        a_m1 = r_m0'*r_m0/(p_m0'*A*p_m0);
        x_m1 = x_m0 + a_m1*p_m0;
        r_m1 = r_m0 - a_m1*A*p_m0;
        b_m1 = r_m1'*r_m1/(r_m0'*r_m0);
        p_m1 = r_m1 + b_m1*p_m0;

        p_m0 = p_m1;
        r_m0 = r_m1;
        x_m0 = x_m1;

        x_stars(:,i) = x_m1;
    end
end