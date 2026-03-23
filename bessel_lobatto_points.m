function [r, x] = bessel_lobatto_points(N, nu, R)
% Computes Bessel–Lobatto points based on zeros of J_nu'(x)
% N  – number of points (including endpoints)
% nu – Bessel order (usually nu = 0)
% R  – radius of domain
%
% r – physical Bessel–Lobatto points in [0, R]
% x – zeros of derivative J_nu'(x) (unscaled)

    if nargin < 3
        R = 1; % Default radius
    end

    % Number of interior points
    M = N - 2;

    % Allocate
    x = zeros(M+1,1);

    % Initial guesses for derivative zeros
    % Use asymptotic approximation for nth zero of J_nu'
    k = (1:M).';
    x0 = (k + nu/2 - 3/4)*pi;   % good first guess

    % Solve Jnu'(x) = 0 using fzero
    for j = 1:M
        f = @(t) besselj(nu-1,t) - (nu/t)*besselj(nu,t);  % J?'(x) relation
        x(j) = fzero(f, x0(j));
    end

    % Append derivative zero used for scaling
    % Use asymptotic approximation for (M+1)th zero, refine
    x_guess = (M+1 + nu/2 - 3/4)*pi;
    f = @(t) besselj(nu-1,t) - (nu/t)*besselj(nu,t);
    x(M+1) = fzero(f, x_guess);

    % Scale to physical domain [0, R]
    r = [0;  R * x(1:M) ./ x(M+1);  R];
end


