%% 3D Cylindrical Magnetorquer Solver using Bessel Functions
% Computes B_r and B_z for a single loop using Bessel function integral

clear; clc;

%% Physical constants
mu0 = 4*pi*1e-7;  % H/m

%% Coil parameters
I = 1.0;          % current [A]
Rcoil = 0.05;     % loop radius [m]
zcoil = 0.0;      % axial position [m]

%% Observation points
rmax = 0.2; Nr = 100;
zmax = 0.3; Nz = 150;
r = linspace(0,rmax,Nr);
z = linspace(-0.05,zmax,Nz);
[R,Z] = meshgrid(r,z);

%% Bessel function integration parameters
Nint = 200;          % number of integration points
kmax = 200;          % maximum Bessel function argument
k = linspace(0,kmax,Nint);  % integration variable

%% Preallocate B fields
Br = zeros(Nz,Nr);
Bz = zeros(Nz,Nr);

%% Compute fields using the analytical formula (see Smythe "Static and Dynamic Electricity")
for i = 1:Nz
    for j = 1:Nr
        r_obs = R(i,j);
        z_obs = Z(i,j) - zcoil;
        
        % B_r and B_z integral using Bessel functions
        integrand_Bz = k .* besselj(0,k*r_obs) .* exp(-k*abs(z_obs)) .* besselj(1,k*Rcoil);
        integrand_Br = k .* besselj(1,k*r_obs) .* exp(-k*abs(z_obs)) .* besselj(1,k*Rcoil);
        
        Bz(i,j) = mu0*I/(2*pi) * trapz(k,integrand_Bz);
        Br(i,j) = mu0*I/(2*pi) * sign(z_obs) * trapz(k,integrand_Br);
    end
end

%% Plot B_z
figure;
contourf(R,Z,Bz,50,'LineColor','none');
colorbar;
xlabel('r [m]');
ylabel('z [m]');
title('Axial Magnetic Field B_z [T]');

%% Plot B_r
figure;
contourf(R,Z,Br,50,'LineColor','none');
colorbar;
xlabel('r [m]');
ylabel('z [m]');
title('Radial Magnetic Field B_r [T]');