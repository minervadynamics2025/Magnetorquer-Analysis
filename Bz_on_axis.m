function [ Bz ] = Bz_on_axis(z, r0, l, J0, mu0)
% Axial magnetic field at r = 0

zp = z + l/2;
zm = z - l/2;

Bz = (mu0 * J0 * r0^2 / 4) * ...
    ( zp ./ (r0^2 + zp.^2).^(3/2) - ...
      zm ./ (r0^2 + zm.^2).^(3/2) );
end

