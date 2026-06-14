function A=alpha_example(alpha,varargin)
% Usage:
%  A=alpha_example(alpha)   or
%  A=alpha_example(alpha,m)
%
% Generates a matrix of size m (which is 20 unless 
% it is specified). alpha is a parameter which should be
% in the interval (0,infty).
%
% The basic QR-method will be have in very different ways
% for large and small alpha
%


if (~(alpha>0))
    alpha
    error('alpha must be positive');
end
if (length(varargin)==0)
    m=10;
else
    m=varargin{1};
end

randn('state',0);
a=2;
d=(a.^(0:(m-1)));
d(end-1)=d(end)*(0.99-1/(5*alpha));
B=randn(m);
A=B\(diag(d)*B);

