function [ out ] = getFieldMetrics( geo, mu_r, I, s )
%GETFIELDMETRICS Core-aware transfer factor and demagnetizing factor.
%
%   out.kappa : <B_z>_rod / B_z(0, z_sensor)  (converts the teslameter reading
%               into the volume-averaged internal flux density)
%   out.Nd    : demagnetizing factor consistent with the volume averages,
%               Nd = (mu0*mu_r*H_app/<B> - 1)/(mu_r - 1)
%   out.Bvol  : volume-averaged B_z in the rod [T]
%   out.Bs    : centreline B_z at the sensor [T]
%   out.Bpeak, out.Bend : centreline B_z at the midplane / rod end face [T]
%   out.z, out.Bax      : centreline profile (z >= 0)
%
%   s (optional) overrides geo.s (sensor distance from the end of the winding).

if nargin < 4, s = geo.s; end
geo.s = s;
mu0 = 4e-7*pi;
[r, z, psi] = solveAxisymmetricField(geo, mu_r, I);

Bax  = 2*psi(2,:)/r(2)^2;                     % centreline field
[~, i0] = min(abs(r - geo.r0));
Bsec = 2*psi(i0,:)/r(i0)^2;                   % cross-section average in the rod
m = z <= geo.L_rod/2 + 1e-12;
Bvol = trapz(z(m), Bsec(m))/(geo.L_rod/2);
Bs = interp1(z, Bax, geo.L_wind/2 + s);
Happ = geo.N*I/geo.L_wind;

out.z = z; out.Bax = Bax;
out.Bvol = Bvol; out.Bs = Bs; out.kappa = Bvol/Bs;
out.Bpeak = Bax(1);
out.Bend = interp1(z, Bax, geo.L_rod/2 + 1e-4);
out.Mmid = interp1(z, Bsec, 0.8*geo.L_rod/2)/Bsec(1);   % M(z)/M(0) at |z| = 0.4 l
out.Mend = interp1(z, Bsec, geo.L_rod/2)/Bsec(1);       % ... at the end face
if mu_r > 1
    out.Nd = (mu0*mu_r*Happ/Bvol - 1)/(mu_r - 1);
else
    out.Nd = NaN;
end
end
