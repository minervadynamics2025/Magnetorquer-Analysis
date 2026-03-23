function [r, z, R, Z] = cylindrical_LGL_Bessel_grid(Nz, Nr, Rmax)
%--------------------------------------------------------------------------
% Generates tensor-product grid:
%   z-direction: Legendre–Gauss–Lobatto (LGL)
%   r-direction: Bessel–Lobatto nodes (zeros of J0' plus endpoints)
%
% INPUTS:
%   Nz   = number of LGL points in z
%   Nr   = number of Lobatto–Bessel points in r
%   Rmax = outer radius
%
% OUTPUTS:
%   r = 1D Lobatto–Bessel nodes in [0, Rmax]
%   z = 1D LGL nodes in [-1, 1]
%   R,Z = 2D meshgrid
%--------------------------------------------------------------------------
    
    % ------------------------------
    % 1. Generate LGL points (z-direction)
    % ------------------------------
    % z = legendre_gauss_lobatto(Nz);
    [z,w] = legendre_gauss_lobatto(Nz);

    % ------------------------------
    % 2. Generate Bessel–Lobatto points (r-direction)
    % ------------------------------
    nu = 0; % Zeroth order Bessel function
    [r, w] = bessel_lobatto_quadrature(Nr, nu, Rmax);
    %r = bessel_lobatto_points(Nr, Rmax);

    % 3. Tensor-product mesh
    [Z, R] = meshgrid(z, r);

end

