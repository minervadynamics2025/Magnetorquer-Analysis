function [ H_f, M ] = getHfieldLangevinModel( N_d, H_0, M_s, a )
%GETHFIELDLANGEVINMODEL Self-consistent magnetization of a rod in an applied field.
%   Solves  M = M_s * L( (H_0 - N_d*M)/a ),  L(x) = coth(x) - 1/x,
%   and returns the effective magnetizing field H_f = H_0 - N_d*M and M.
%
% N_d : Demagnetization factor
% H_0 : Applied magnetizing field H_app [A/m]
% M_s : Saturation magnetization [A/m]
% a   : Shape parameter [A/m]
%
% The root is bracketed in (0, min(H_0/N_d, M_s)): H_eff must stay positive and
% M cannot exceed M_s, so fzero always converges (the previous version used an
% unbracketed initial guess, which could converge to a non-physical root).

L = @(x) coth(x) - 1./x;
f = @(x) x - M_s * L((H_0 - N_d*x)/a);

upper = min(H_0/N_d, M_s) * (1 - 1e-12);
M = fzero(f, [1e-9, upper]);
H_f = H_0 - N_d*M;
end
