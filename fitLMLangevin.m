function [ params ] = fitLMLangevin( B, H_0, N_d, params )
%% LM fitting of Langevin model to magnetorquer B-field measurements
% N_d : Demagnetization factor
% H_0 : Applied magnetazing field
% B   : Measured Magnetic Field
% params: Initial guess from Simulated Annealing

M_s = params(1);         % Approximate saturation magnetization
a_0 = params(2);         % Shape Parameter
mu0 = 4*pi*1e-7;         % Permeability

%% Langevin model function
syms Ms H a

% Langevin function analytically
B_sym = mu0 * (H + Ms * (coth(H/a) - a/H));        % total B-field

% Convert symbolic function to MATLAB function handle for numeric evaluation
B_func = matlabFunction(B_sym, 'Vars', [Ms, a, H]);

% Symbolic derivatives (Jacobian)
dB_dMs = matlabFunction(diff(B_sym, Ms), 'Vars', [Ms, a, H]);
dB_dAlpha = matlabFunction(diff(B_sym, a), 'Vars', [Ms, a, H]);

%% Initial guess for [Ms, alpha]
params0 = [M_s, a_0];  % Ms in Tesla, alpha unitless

%% Levenberg-Marquardt settings
maxIter = 1000;
lambda = 0.01;
tol = 1e-8;

B_calc = zeros(length(H_0),1); % Model Magnetic Field
H_f = zeros(length(H_0),1);    % Effective Magnetizing Field
params = params0;
for k = 1:maxIter
    
    % Jacobian using symbolic derivatives
    J = zeros(length(H_0), 2);
    % Here we compute corresponding fields
    for i=1:length(H_0)
        H_f(i) = getHfieldLangevinModel( N_d, H_0(i), params(1), params(2) );
        B_calc(i) = B_func(params(1), params(2), H_f(i));
        J(i,1) = dB_dMs(params(1), params(2), H_f(i));
        J(i,2) = dB_dAlpha(params(1), params(2), H_f(i));
    end     
    % Compute residuals
    r = B - B_calc;
    
    % LM update
    delta_p = (J.'*J + lambda*eye(2)) \ (J.'*r);
    %delta_p = (J.'*J) \ (J.'*r);
    params_new = params + delta_p';
    
    % Here we compute corresponding magnetic field with new parameters
    for i=1:length(H_0)
        B_calc(i) = B_func(params_new(1), params_new(2), H_f(i));
    end   
    
    %Check residual
    r_new = B - B_calc;
    
    if norm(r_new) < norm(r)
        params = params_new;
        lambda = lambda / 10;
    else
        lambda = lambda * 10;
    end
    norm(delta_p)
    % Convergence
    if norm(delta_p) < tol
        break;
    end
end


end