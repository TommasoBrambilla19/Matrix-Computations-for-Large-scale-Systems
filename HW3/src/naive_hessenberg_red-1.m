function H=naive_hessenberg_red(A);
% A naive inefficient way of carrying out 
% the hessenberg reduction 
%

n=length(A);

for k=1:n-2
    x=A(k+1:n,k);
    alpha = -sign(x(1))*norm(x);
    z=x-alpha*[1;zeros(length(x)-1,1)];
    u=z/norm(z);

    P1=eye(n-k)-2*u*u';
    P0=eye(k);
    P=[P0,zeros(k,n-k); zeros(n-k,k),P1];  % Equation (2.5) in
                                           % lecture notes

    PA=P*A;  
    % PA should have zeros below the second off-diagonal in column
    % k as described on page 6 in lecture notes.
    A=PA*P';
    % A should have the same zero-structure as PA 

end
H=A; % should be a hessenberg matrix with same eigenvalues as input A
