function [ mu_c, rms ] = constantPermeabilityFit( B, H_0, N_d )
%CONSTANTPERMEABILITYFIT Best constant mu_r for <B> = mu0*mu_r*H_0/(1 + N_d(mu_r - 1)).
% Reference model: in the demagnetization-limited regime every mu_r >> 1/N_d
% gives nearly the same line mu0*H_0/N_d, so this fit is the benchmark the
% Langevin model has to beat.
mu0 = 4e-7*pi;
model = @(m) mu0*m*H_0(:)./(1 + N_d*(m - 1));
cost = @(lm) sum(((model(exp(lm)) - B(:))./B(:)).^2);
lm = fminsearch(cost, log(3000));
mu_c = exp(lm);
rms = sqrt(cost(lm)/numel(B));
end
