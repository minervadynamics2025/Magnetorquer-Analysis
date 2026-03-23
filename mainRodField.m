% Here we solve the magnetic potential equation with spectral methods
clear; clc;

a = 1;        % domain radius
m = 0;        % axisymmetric
N = 10;       % number of Bessel modes
Nr = 500;     % resolution

fhandle = @(r) exp(-5*r);   % example function

[r, f_approx, coeffs, alpha] = bessel_series_zeros(fhandle, a, m, N, Nr);

plot(r, fhandle(r), 'k', 'LineWidth', 2); hold on;
plot(r, f_approx, '--r', 'LineWidth', 1.5);
legend('Original f(r)', 'Bessel Series Approx.');
xlabel('r'); ylabel('f(r)');
title('Fourier–Bessel Series Representation');

% Compute Legendre-Gauss-Lobatto points
N = 6;
[x,w] = legendre_gauss_lobatto(N);

% Display
fprintf('LGL nodes (N=%d):\n',N); disp(x')
fprintf('LGL weights (sum should be 2):\n'); disp(w')
fprintf('Sum of weights = %g\n', sum(w));

% Compute Bessel–Lobatto Points and Weights

N = 12;      % total points
nu = 0;      % J0
R = 1.0;

[r, w] = bessel_lobatto_quadrature(N, nu, R);

% plot nodes and weights
disp(table((0:N-1).', r, w, 'VariableNames', {'idx','r','w'}));
figure; plot(r, w, 'o-'); xlabel('r'); ylabel('w'); title('Bessel-Lobatto weights');

% Create a grid with LGL points in z direction
% and Bessel-Lobatto points in r direction
Nz  = 3;      % LGL points in z
Nr  = 5;      % Bessel-Lobatto points in r
Rmax = 1.0;    % max radius

[r, z, R, Z] = cylindrical_LGL_Bessel_grid(Nz, Nr, Rmax);

figure; plot(z, zeros(size(z)), 'o'); title('LGL points');
figure; plot(r, zeros(size(r)), 'o'); title('Lobatto-Bessel points');

figure; scatter(Z(:), R(:), 'x'); xlabel('z'); ylabel('r');
title('Grid Points: LGL (z) × Bessel-Lobatto (r)');

%% Representing a Bessel Function Using Lagrange Interpolants
nu = 0;
N = 20;
% Lobatto-Bessel or simply scaled Bessel-derivative zeros; here use linspace for demo
r_nodes = linspace(0,10,N).';
r_eval = linspace(0,10,1000).';
f_interp = bessel_lagrange_barycentric(nu, r_nodes, r_eval);
f_true = besselj(nu, r_eval);

figure;
plot(r_eval, f_true, 'k', 'LineWidth', 1.5); hold on;
plot(r_eval, f_interp, '--r', 'LineWidth', 1);
plot(r_nodes, besselj(nu, r_nodes), 'ob', 'MarkerFaceColor','b');
legend('true J_\nu', 'barycentric interpolant', 'nodes');
xlabel('r'); ylabel('J_\nu(r)');
title('Barycentric Lagrange Interpolation of Bessel function');

