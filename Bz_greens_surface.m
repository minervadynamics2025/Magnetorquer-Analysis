function Bz = Bz_greens_surface(r, z, r0, l, J0, mu0)
% Computes Bz using Green's function on surface z = l/2
%
% Inputs:
%   r         - observation radial position (scalar or array)
%   z         - observation axial position
%   r0        - radius of surface
%   l         - solenoid length
%   sigma_fun - function handle: sigma(r')
%   J0  - current density
%   mu0 - permeability (e.g., 4*pi*1e-7)
%
% Output:
%   Bz        - magnetic field

Nr = 20;
% Integration grid
rp = linspace(0, r0, Nr);
dr = rp(2) - rp(1);

% Surface location
zp = l/2;

Bs = zeros(1,length(rp));  % Here we initialize Bz vector
% Evaluate Magnetic Field at Surface
for i=1:length(rp)
    Bs(i) = Bz_analytical(rp(i), l/2, r0, l, J0, mu0);
end

% Plot
figure;
plot(rp, Bs, 'LineWidth', 2);
grid on;

% Distance
R = sqrt((rp-r).^2 + z^2);

% Kernel
kernel = z * dr ./ (2.* R.^3.);

% Integrand
integrand = Bs .* rp .* kernel;

% Integrate
Bz = trapz(rp, integrand);


end