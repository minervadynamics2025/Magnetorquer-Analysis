function [ H_f ] = getHfieldLangevinModel( N_d, H_0 )
%This function calculates the root of hyperbolic tangent model
% N_d : Demagnetization factor
% H_0 : Applied magnetazing field
M_s = 4.61 * 1e3; % Approximate saturation magnetization
a = 600.;         % Shape Parameter

% Hyperbolic tangent model
f = @(x) x -  M_s * (coth((H_0 - N_d * x)/a) - a/(H_0 - N_d * x));

% Initial guess
x0 = H_0;

% Find root
M = fzero(f, x0);

% Effective magnetizing field
H_f = H_0 - N_d * M;

end

