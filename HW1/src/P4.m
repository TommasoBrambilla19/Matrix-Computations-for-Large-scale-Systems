load('Bwedge.mat');

%% Part A
% Find the three eigenvalues that are furthest from each other
distances = abs(B_eigvals - B_eigvals.');
norms = vecnorm(distances);
[~, idx] = maxk(norms,3);
% Plot the eigenvalues
scatter(real(B_eigvals), imag(B_eigvals), 'filled');
hold on;
plot(real(B_eigvals(idx(1))), imag(B_eigvals(idx(1))), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(B_eigvals(idx(2))), imag(B_eigvals(idx(2))), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(B_eigvals(idx(3))), imag(B_eigvals(idx(3))), 'gx', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;

%% Part B
%%% First point(red cross)
R_eigvals = B_eigvals(setdiff(1:length(B_eigvals), idx));
[r1, c1] = neat_plot(R_eigvals, 'r', [-10, 0]);
e1 = [real(B_eigvals(idx(1))) imag(B_eigvals(idx(1)))];
conv_fact_1 = r1/norm(e1-c1)

%%% Second point(black cross)
K_eigvals = B_eigvals(setdiff(1:length(B_eigvals), idx));
[r2, c2] = neat_plot(K_eigvals, 'k', [-10, -30]);
e2 = [real(B_eigvals(idx(2))) imag(B_eigvals(idx(2)))];
conv_fact_2 = r2/norm(e2-c2)

%%% Third point(green cross)
G_eigvals = B_eigvals(setdiff(1:length(B_eigvals), idx));
[r3, c3] = neat_plot(G_eigvals, 'g', [-10, 30]);
e3 = [real(B_eigvals(idx(3))) imag(B_eigvals(idx(3)))];
conv_fact_3 = r3/norm(e3-c3)
hold off

%% Part C
figure; % Create a new figure for Part C
b = ones(length(B(:, 1)), 1);
colors = lines(7); % Generate 7 unique colors
m_values = [2, 4, 8, 10, 20, 30, 40];
initial_size = 80; % Initial size of the markers

err_vals_1 = zeros(1, length(m_values));
err_vals_2 = zeros(1, length(m_values));
err_vals_3 = zeros(1, length(m_values));

for i = 1:length(m_values)
    m = m_values(i);
    [Q, H] = arnoldi(B, b, m);
    ansArn = eig(H(1:m, 1:m));
    scatter(real(ansArn), imag(ansArn), initial_size - (i-1)*10, colors(i, :), 'filled'); % Use unique color for each m and decrease size
    hold on;

    err_vals_1(i) = abs(min(ansArn-B_eigvals(idx(1))));
    err_vals_2(i) = abs(min(ansArn-B_eigvals(idx(2))));
    err_vals_3(i) = abs(min(ansArn-B_eigvals(idx(3))));
end
hold off;

figure;
semilogy(m_values, err_vals_1, 'r', 'DisplayName', 'Red cross');
hold on;
semilogy(m_values, conv_fact_1.^m_values, "k--")
figure;
semilogy(m_values, err_vals_2, 'k', 'DisplayName', 'Black cross');
hold on;
semilogy(m_values, conv_fact_2.^m_values, "k--")
figure;
semilogy(m_values, err_vals_3, 'g', 'DisplayName', 'Green cross');
hold on;
semilogy(m_values, conv_fact_3.^m_values, "k--")
hold off

%% Part D
ref_eigval = B_eigvals(12);
b = ones(length(B(:, 1)), 1);
shifts = [-10, -7+2i, -9.8+1.5i];
m_values = [10, 20, 30];

% Initialize table to store results
results = table('Size', [length(shifts)*length(m_values), 3], ...
                'VariableTypes', {'double', 'double', 'double'}, ...
                'VariableNames', {'Shift', 'm', 'Error'});

row = 1;
for shift = shifts
    for m = m_values
        C = inv(B - shift * eye(size(B)));
        [Q, H] = arnoldi(C, b, m);
        lambdas = shift + 1 ./ eig(H(1:m, 1:m));
        [~, i] = min(abs(lambdas - ref_eigval));
        closest_eigval = lambdas(i);
        lambda_error = abs(closest_eigval - ref_eigval);
        
        % Store results in the table
        results.Shift(row) = shift;
        results.m(row) = m;
        results.Error(row) = lambda_error;
        row = row + 1;
    end
end

% Display the results table
disp(results);
hold off;
%% Copied code from stackoverflow https://stackoverflow.com/questions/17377474/defining-a-minimal-bounding-circle-in-matlab
function r = radiusFromPoint(P, pointsX, pointsY)
    px = P(1);
    py = P(2);
    distanceSquared = (pointsX - px).^2 + (pointsY - py).^2;
    r = sqrt(max(distanceSquared));
end

function [r, c]=neat_plot(C_eigvals, color, P0)
    pointsX = real(C_eigvals); % Copied from stackoverflow
    pointsY = imag(C_eigvals);
    [P, radiusMin] = fminsearch(@(P) radiusFromPoint(P, pointsX, pointsY), P0);

    theta = linspace(0, 2*pi, 100);
    x_circle = P(1) + radiusMin * cos(theta);
    y_circle = P(2) + radiusMin * sin(theta);
    plot(x_circle, y_circle, [color '--']);
    plot(P(1), P(2), [color '.']);

    r = radiusMin;
    c = P;
end
