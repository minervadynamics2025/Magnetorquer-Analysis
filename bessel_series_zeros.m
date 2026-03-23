function [r, f_approx, coeffs, alpha] = bessel_series_zeros(fhandle, a, m, N, Nr)
% bessel_series_zeros
% Computes Fourier–Bessel expansion of order m using the
% first N zeros of J_m.
%
% Inputs:
%   fhandle : function handle f(r)
%   a       : radius domain [0, a]
%   m       : Bessel order (usually 0 for axisymmetric)
%   N       : number of Bessel modes
%   Nr      : number of radial sample points
%
% Outputs:
%   r        : radial grid
%   f_approx : reconstructed series approximation
%   coeffs   : Bessel coefficients c_n
%   alpha    : zeros of J_m

    % --- radial grid ---
    r = linspace(0, a, Nr);
    dr = r(2) - r(1);
    w = r;                   % weight r dr for integral

    % --- compute zeros of J_m ---
    alpha = besselzero(m, N);

    % --- evaluate target function ---
    f = fhandle(r);

    % --- compute coefficients ---
    coeffs = zeros(1, N);

    for n = 1:N
        % mode shape
        Jmn = besselj(m, alpha(n)*r/a);

        % numerator: integral f(r) * r * J_m(...)
        num = sum(f .* w .* Jmn) * dr;

        % denominator from orthogonality
        denom = (a^2/2) * (besselj(m+1, alpha(n)))^2;

        coeffs(n) = num / denom;
    end

    % --- reconstruct function ---
    f_approx = zeros(size(r));
    for n = 1:N
        f_approx = f_approx + coeffs(n) * besselj(m, alpha(n)*r/a);
    end
end


