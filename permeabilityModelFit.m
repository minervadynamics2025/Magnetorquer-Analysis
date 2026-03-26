function [ mu_fit] = permeabilityModelFit( mu_r, H_f )
% this function fits relative permeability 
% mu_r : Relative Permeability calculated from measurements 
% H_f  : Effective Magnetizing Field

model = @(p,H) p(1)*exp(-H/p(2));  % p = [a,b,c]

p0 = [3000, 500];   % initial guess

lb = [0, 0];         % lower bounds
ub = [1e4, 1e5];   % upper bounds

p_fit = lsqcurvefit(model, p0, H_f, mu_r, lb, ub);

% Evaluate fitted curve
mu_fit = model(p_fit,H_f);

end

