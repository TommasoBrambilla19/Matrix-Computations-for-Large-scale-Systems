clc
clear all

rand('seed', 0);
nn = 1200;
m = 100;
A = gallery('wathen', nn, nn);
b = ones(length(A(:,1)),1);

for m=[100, 50, 20, 10, 5]
    [Q, H] = arnoldi(A, b, m);
    should_be_zero1=norm(Q*H-A*Q(:,1:m))
    should_be_zero2=norm(Q'*Q-eye(m+1))
end

% Modified GS: nn=1200, time = 55.848073s, m = 100, Indicator = 2.0351e-12
% Modified GS: nn=1200, time = 17.115119s, m = 50, Indicator = 4.5778e-13
% Modified GS: nn=1200, time = 4.166520s, m = 20, Indicator = 5.7906e-14
% Modified GS: nn=1200, time = 1.683719s, m = 10, Indicator = 3.9879e-13
% Modified GS: nn=1200, time = 0.708851s, m = 5, Indicator = 3.9879e-13

% Triple GS: nn=750, time = 58.430678s, m = 100, Indicator = 4.3691e-13
% Triple GS: nn=750, time = 16.081235s, m = 50, Indicator = 1.2598e-13
% Triple GS: nn=750, time = 3.244074s, m = 20, Indicator = 2.4883e-14
% Triple GS: nn=750, time = 1.057311s, m = 10, Indicator = 1.003e-13
% Triple GS: nn=750, time = 0.415080s, m = 5, Indicator = 1.000e-13

% Double GS: nn=900, time = 59.777176s, m = 100, Indicator = 1.5851e-12
% Double GS: nn=900, time = 16.851038s, m = 50, Indicator = 4.9644e-13
% Double GS: nn=900, time = 3.551333s, m = 20, Indicator = 8.9858e-14
% Double GS: nn=900, time = 1.252394s, m = 10, Indicator = 1.0202e-13
% Double GS: nn=900, time = 0.486280s, m = 5, Indicator = 9.7917e-14

%Single GS: nn=1200, time = 56.246460, m = 100, Indicator = 3.8425e-12
%Single GS: nn=1200, time = 17.380391, m = 50, Indicator = 1.0349e-12
%Single GS: nn=1200, time = 4.272772, m = 20, Indicator = 1.4068e-13
%Single GS: nn=1200, time = 1.657926, m = 10, Indicator = 3.9879e-13
%Single GS: nn=1200, time = 0.731169, m = 5, Indicator = 3.9879e-13