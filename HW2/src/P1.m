%% B
alpha = [1, 5, 10, 100];
n=100;
rand('state',5);
b = rand(n,1);
A = sprand(n,n,0.5);
m = 100;

figure
conv = zeros(2, length(alpha));

for i = 1:length(alpha)
    A1 = A + alpha(i)*speye(n); 
    [C_evec, C_eval] = eig(full(A1));

    C_eigs = diag(C_eval);
    K_eigs = C_eigs(setdiff(1:length(C_eigs), 1));

    
    subplot(2, 2, i)
    scatter(real(C_eigs), imag(C_eigs), 'kx')
    hold on
    [r, c] = neat_plot(K_eigs, 'r', [alpha, 0]);
    c = complex(c(1), c(2));

    title(['\alpha = ', num2str(alpha(i))])
    xlabel('Real Part')
    ylabel('Imaginary Part')
    axis equal
    grid on
    legend('Eigvals','Localization disk', 'Center point', 'Location', 'best');
    hold off
    conv(:, i) = [norm(C_evec)*norm(inv(C_evec))*(norm(C_eigs(1)-c) + r)/norm(C_eigs(1)); r/norm(c)];
end

%% A
alpha = [1, 5, 10, 100];
n=100;
rand('state',5);
b = rand(n,1);
A = sprand(n,n,0.5);
m = 100;

figure
for i = 1:length(alpha)   
    A1 = A + alpha(i)*speye(n); 
    A1=A1/norm(A1,1);

    [Q, H, x_stars, z_stars, ~]= gmres_arnoldi(A1,b,m);
    

    relative_error = vecnorm(x_stars-x_stars(:,end)) / norm(x_stars(:,end));
    residual_norm = vecnorm(A1*x_stars-b) / norm(b);

    subplot(2, 2, i)
    semilogy(relative_error, 'r-')
    hold on
    semilogy(residual_norm, 'k--')
    semilogy(conv(1,i)*conv(2, i).^(0:m-1), 'b--')
    hold off
    legend('|Ax-b|/|b|', '|x-x_{end}|/|x_{end}|', 'theory prediction', 'Location', 'best')
    title(['\alpha = ', num2str(alpha(i))])
    ylim([1e-17, 1e5])
end

%% C 
results = [];
results2 = [];
for n = [200 500 1000]
    b = rand(n,1);
    A = sprand(n,n,0.5);
    for a = [1 100]
        A1 = A + a*speye(n);
        A1 = A1/norm(A1,1);

        tic;
        x_stars_2 = A1\b;
        elapsedTime2 = toc;

        resnorm2 = norm(A1*x_stars_2-b);

        results2 = [results2; a, n, elapsedTime2, resnorm2];

        for m = [5 10 20 50 100]
            tic;
            [Q, H, x_stars, z_stars, ~] = gmres_arnoldi(A1, b, m);
            elapsedTime = toc;

            resnorm = norm(A1*x_stars(:,end)-b);

            results = [results; a, m, n, elapsedTime, resnorm];


        end
    end
end

for i = 1:size(results, 1)
    fprintf('alpha = %d, m = %d, n = %d, Time = %.4f, Resnorm = %d. \n', results(i, 1), results(i, 2), results(i, 3), results(i, 4), results(i, 5));
end

disp('\n')
for i = 1:size(results2, 1)
    fprintf('alpha = %d, n = %d, Time = %.4f, Resnorm = %d. \n', results2(i, 1), results2(i, 2), results2(i, 3), results2(i, 4));
end

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
    plot(x_circle, y_circle, [color '-'], "DisplayName", join(["r = ", num2str(radiusMin)]));
    plot(P(1), P(2), [color '.'], "DisplayName", join(["c = [",num2str(P(1)), num2str(P(2)) "]"]));

    r = radiusMin;
    c = P;
end
