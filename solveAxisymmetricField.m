function [ r, z, psi ] = solveAxisymmetricField( geo, mu_r, I )
%SOLVEAXISYMMETRICFIELD Axisymmetric magnetostatic field of a coil wound on a
% ferromagnetic rod, with the core included (finite-volume method).
%
%   Solves  d/dr( 1/(mu r) d(psi)/dr ) + d/dz( 1/(mu r) d(psi)/dz ) = -J_theta
%   for the flux function psi = r*A_theta on the half domain z >= 0
%   (symmetry plane z = 0, psi = 0 on the axis and on the far boundary).
%   B_z = (1/r) d(psi)/dr, and the flux through a disc of radius r is 2*pi*psi.
%
%   geo  : struct with fields
%          N      - number of turns
%          L_wind - winding length [m]
%          L_rod  - total rod length [m] (>= L_wind)
%          r0     - rod radius [m]
%          t_coil - radial build of the winding [m]
%          r_in   - (optional) inner radius of the winding [m], default r0
%          s      - sensor distance from the end of the winding [m]
%   mu_r : relative permeability of the core (1 = air core)
%   I    : coil current [A]
%
%   Returns grid vectors r, z and psi(nr, nz).

mu0 = 4e-7*pi;

% ---------------- graded grid
if ~isfield(geo, 'r_in'), geo.r_in = geo.r0; end   % inner radius of the winding
r = linspace(0, geo.r0, 41);
if geo.r_in > geo.r0
    rg = linspace(geo.r0, geo.r_in, 9); r = [r, rg(2:end)];
end
rc = linspace(geo.r_in, geo.r_in + geo.t_coil, 21); r = [r, rc(2:end)];
dr = r(end) - r(end-1);
while r(end) < 0.6
    dr = dr*1.08; r(end+1) = r(end) + dr; %#ok<AGROW>
end
zend = max(geo.L_rod/2, geo.L_wind/2 + geo.s + 0.002) + 0.012;
z = linspace(0, geo.L_wind/2, 301);
nfine = floor((zend - geo.L_wind/2)/1.25e-4) + 1;
z2 = linspace(geo.L_wind/2, zend, nfine);
z = [z, z2(2:end)];
dz = z(end) - z(end-1);
while z(end) < 0.7
    dz = dz*1.08; z(end+1) = z(end) + dz; %#ok<AGROW>
end
nr = numel(r); nz = numel(z);
drs = diff(r); dzs = diff(z);
rcc = 0.5*(r(1:end-1) + r(2:end));
zcc = 0.5*(z(1:end-1) + z(2:end));
[RC, ZC] = ndgrid(rcc, zcc);

% ---------------- cell properties
nu = ones(size(RC));                          % 1/mu_r
nu(RC < geo.r0 & ZC < geo.L_rod/2) = 1/mu_r;
J = zeros(size(RC));
J(RC > geo.r_in & RC < geo.r_in + geo.t_coil & ZC < geo.L_wind/2) = ...
    geo.N*I/(geo.L_wind*geo.t_coil);

% ---------------- assembly (node-centred control volumes)
idx = @(i,j) (j-1)*nr + i;                    % column-major numbering
nmax = 5*nr*nz;
rows = zeros(nmax,1); cols = rows; vals = rows; k = 0;
b = zeros(nr*nz,1);
for j = 1:nz
    for i = 1:nr
        p = idx(i,j);
        if i == 1 || i == nr || j == nz       % Dirichlet psi = 0
            k = k+1; rows(k) = p; cols(k) = p; vals(k) = 1; continue
        end
        jm = max(j-1, 1);
        if j > 1, dzm = dzs(j-1); else, dzm = 0; end
        dzp = dzs(j);
        hm = dzm/2; hp = dzp/2;
        diagv = 0;
        % east / west neighbours
        a = (nu(i,jm)*hm + nu(i,j)*hp)/(rcc(i)*drs(i));
        k = k+1; rows(k) = p; cols(k) = idx(i+1,j); vals(k) = a; diagv = diagv - a;
        a = (nu(i-1,jm)*hm + nu(i-1,j)*hp)/(rcc(i-1)*drs(i-1));
        k = k+1; rows(k) = p; cols(k) = idx(i-1,j); vals(k) = a; diagv = diagv - a;
        % north / south neighbours
        ra = r(i) - drs(i-1)/4; rb = r(i) + drs(i)/4;
        g = nu(i-1,j)*(drs(i-1)/2)/ra + nu(i,j)*(drs(i)/2)/rb;
        k = k+1; rows(k) = p; cols(k) = idx(i,j+1); vals(k) = g/dzp; diagv = diagv - g/dzp;
        if j > 1
            g = nu(i-1,j-1)*(drs(i-1)/2)/ra + nu(i,j-1)*(drs(i)/2)/rb;
            k = k+1; rows(k) = p; cols(k) = idx(i,j-1); vals(k) = g/dzm; diagv = diagv - g/dzm;
        end
        k = k+1; rows(k) = p; cols(k) = p; vals(k) = diagv;
        src = (J(i-1,j)*drs(i-1)/2 + J(i,j)*drs(i)/2)*hp;
        if j > 1
            src = src + (J(i-1,j-1)*drs(i-1)/2 + J(i,j-1)*drs(i)/2)*hm;
        end
        b(p) = -mu0*src;
    end
end
A = sparse(rows(1:k), cols(1:k), vals(1:k), nr*nz, nr*nz);
psi = reshape(A\b, nr, nz);
end
