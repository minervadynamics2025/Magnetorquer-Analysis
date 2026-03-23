%% Magnetic Dipole Moment with Saturation for Hiperco 50A Rod
clear; clc;

%% Parameters
mu0 = 4*pi*1e-7;          % Vacuum permeability [H/m]
mu_r = 3000;              % Relative permeability of Hiperco 50A (linear region)
B_sat = 2.35;             % Saturation flux density [T]
l_mag = 0.125;             % Core length [m]
I = 0.06;                  % Coil current [A]
d_copper = 0.0003;       % Diameter of copper wire
r_copper = 1.72e-8;      % Resistivity of copper

% Range of coil parameters
N_array = 400:5:800;         % Number of turns
d_array = 0.001:0.0005:0.005; % Core diameter [m]

%% Preallocate results
m_results = zeros(length(N_array), length(d_array));
R_array = zeros(length(N_array), length(d_array));
V_array = zeros(length(N_array), length(d_array));


%% Loop over turns and diameter
for i = 1:length(N_array)
    N = N_array(i);
    for j = 1:length(d_array)
        d = d_array(j);
        A = pi*(d/2)^2;           % Cross-sectional area [m^2]
        V = A * l_mag;            % Core volume [m^3]
        H = N * I / l_mag;        % Magnetizing field [A/m]

        % Linear flux density
        B_linear = mu_r * mu0 * H;

        % Limit flux to saturation
        B_effective = min(B_linear, B_sat);

        % Dipole moment
        m_results(i,j) = B_effective * V / mu0;
        
        % Resistance of Wire
        R_array(i,j) = (4 * r_copper * N * d) / (d_copper^2);   
        
        % Voltage
        V_array(i,j) = I * R_array(i,j);         
        
    end
end

%% Display results
disp('Magnetic dipole moment (A·m^2) with saturation:');
disp(m_results);

disp('Voltage (V):');
disp(V_array);

%% Surface plot
figure;
[X,Y] = meshgrid(d_array,N_array);
surf(X,Y,m_results)
xlabel('Core Diameter [m]')
ylabel('Number of Turns')
zlabel('Magnetic Dipole Moment [A·m^2]')
title('Dipole Moment with Saturation for Hiperco 50A')
colorbar


% Sinle Dipole Moment Calculation
d = 0.01;                 % Diameter of the rod
N = 2800;                 % Number of turns
A = pi*(d/2)^2;           % Cross-sectional area [m^2]
V = A * l_mag;            % Core volume [m^3]
H = N * I / l_mag;        % Magnetizing field [A/m]
R = (4 * r_copper * N * d) / (d_copper^2);   % Resistance of Wire
V_supply = I * R;         % Voltage

% Dimagnetization factor
N_d = (4. * log(2.*l_mag/d) -1)/((2.*l_mag/d)^2 - 4. * log(2.*l_mag/d)); 
% Magnetization 
M = ((mu_r-1.) * H)/(1+N_d*(mu_r-1)); 

% Dipole moment
m_dipole = M * A * l_mag;
