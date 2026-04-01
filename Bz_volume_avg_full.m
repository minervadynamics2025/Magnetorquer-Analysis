function [Bz_avg] = Bz_volume_avg_full(r0, l, J0, mu0)
% Computes the volume-averaged Bz using the full expression
Nr = 25;
Nz = 100;
% Grid
r = linspace(0, r0, Nr);
z = linspace(0, l, Nz);
[RR, ZZ] = meshgrid(r, z);

% z+ and z-
zp = ZZ + l/2;
zm = ZZ - l/2;

% F+ and F-
Fp = zp ./ sqrt(RR.^2 + r0^2 + zp.^2);
Fm = zm ./ sqrt(RR.^2 + r0^2 + zm.^2);

% First term
T1 = (RR .* (RR.^2 + 3*r0^2)) ./ (RR.^2 + r0^2).^2 .* (Fp - Fm);

% Second term
T2 = (RR.^2 ./ (RR.^2 + r0^2)) .* ...
     ( zp ./ (RR.^2 + r0^2 + zp.^2).^(3/2) - ...
       zm ./ (RR.^2 + r0^2 + zm.^2).^(3/2) );

% Full Bz
Bz = (mu0 * J0 * r0^2 / 4) * (T1 - T2);

% Multiply by r for cylindrical coordinates
Bz_rweighted = Bz .* RR;

% Integrate over r first (trapz handles uneven spacing well)
Bz_rint = trapz(r, Bz_rweighted, 2); % integrate along r (columns)

% Integrate over z
Bz_volume = trapz(z, Bz_rint); 

% Divide by total volume
V = pi * r0^2 * l;
Bz_avg = (2 * pi / V) * Bz_volume; % factor 2*pi comes from theta integration

end

