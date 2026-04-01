function [ Bz ] = Bz_analytical(r, z, r0, l, J0, mu0)
% Bz_analytical computes axial magnetic field Bz(r,z)
% for a finite solenoid using analytical expression
%
% Inputs:
%   r   - radial position (scalar, vector, or matrix)
%   z   - axial position (same size as r or scalar)
%   r0  - solenoid radius
%   l   - solenoid length
%   J0  - current density
%   mu0 - permeability (e.g., 4*pi*1e-7)
%
% Output:
%   Bz  - axial magnetic field

% Ensure element-wise operations
r = r + 0; 
z = z + 0;

% z+ and z-
zp = z + l/2;
zm = z - l/2;

% Common term
R2 = r.^2 + r0^2;

% F+ and F-
Fp = zp ./ sqrt(R2 + zp.^2);
Fm = zm ./ sqrt(R2 + zm.^2);

% First term
term1 = (r .* (r.^2 + 3*r0^2)) ./ (R2.^2) .* (Fp - Fm);

% Second term
term2 = (r.^2 ./ R2) .* ...
    ( zp ./ (R2 + zp.^2).^(3/2) - ...
      zm ./ (R2 + zm.^2).^(3/2) );

% Final expression
Bz = (mu0 * J0 * r0^2 / 4) .* (term1 - term2);

end

