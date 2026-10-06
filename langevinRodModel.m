function [ B, M, H_f ] = langevinRodModel( p, H_0, N_d )
%LANGEVINRODMODEL Volume-averaged flux density of a rod with Langevin magnetization.
%   p   = [M_s, a]
%   B   = mu0*(H_eff + M) = mu0*(H_0 + (1 - N_d)*M)
%   M   : self-consistent magnetization, H_f : effective field (both column vectors)

mu0 = 4e-7*pi;
n = numel(H_0);
M = zeros(n,1); H_f = zeros(n,1);
for i = 1:n
    [H_f(i), M(i)] = getHfieldLangevinModel(N_d, H_0(i), p(1), p(2));
end
B = mu0*(H_f + M);
end
