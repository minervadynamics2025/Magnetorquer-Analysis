function [ Bz ] = Bz_on_axis(z, r0, l, K0, mu0)
%BZ_ON_AXIS On-axis field of a finite AIR-CORE solenoid (thin current sheet).
%   Bz(0,z) = mu0*K0/2 * [ z+/sqrt(r0^2+z+^2) - z-/sqrt(r0^2+z-^2) ],  z+- = z +- l/2
%   K0 = N*I/l is the surface current density [A/m].
%
% NOTE: this is the coil-only field. With a ferromagnetic core the field at an
% external sensor is ~99% due to the core magnetization and has a different
% spatial distribution, so ratios of this field must NOT be used to convert a
% sensor reading into the internal field; use getFieldMetrics.m instead.
% (The previous version used the kernel z/(r0^2+z^2)^(3/2), which is the field
% of a single loop, not of a solenoid.)

zp = z + l/2;
zm = z - l/2;
Bz = (mu0*K0/2) * ( zp./sqrt(r0^2 + zp.^2) - zm./sqrt(r0^2 + zm.^2) );
end
