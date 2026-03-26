function [ H ] = getHFieldHTangentModel( N_d, H_0 )
%This function calculates the root of hyperbolic tangent model
% N_d : Demagnetization factor
% H_0 : Applied magnetazing field
M_s = 3.61 * 1e3; % Approximate saturation magnetization
a = 400.;         % Shape Parameter

% Hyperbolic tangent model
f = @(x) x - H_0 + N_d * M_s * tanh(x/a);

% Initial guess
x0 = H_0;

% Find root
H = fzero(f, x0);


end

