function [ params, info ] = fitLMLangevin( B, H_0, N_d, params, M_s_fixed )
%FITLMLANGEVIN Levenberg-Marquardt refinement of the Langevin parameters.
% B         : volume-averaged flux density converted from the measurements [T]
% H_0       : applied magnetizing field H_app [A/m]
% N_d       : demagnetization factor
% params    : initial guess [M_s, a] (from simulated annealing)
% M_s_fixed : (optional) keep M_s fixed at this value and fit only a
%
% info.cov  : parameter covariance  sigma^2 (G'G)^-1, sigma^2 = S/(n - n_p)
% info.sd   : standard deviations, info.corr : correlation(M_s, a) (free fit)
% info.rms  : RMS relative residual
%
% Changes from the previous version
%  * the Jacobian is computed by finite differences of the FULL model, i.e.
%    including the dependence of H_eff on (M_s, a) through the demagnetizing
%    field (the symbolic Jacobian held H_eff fixed, which is incorrect)
%  * relative residuals; optional fixed M_s; covariance output

if nargin < 5, M_s_fixed = []; end
fixMs = ~isempty(M_s_fixed);
if fixMs, params(1) = M_s_fixed; end
free = [~fixMs, true];

resid = @(p) (langevinRodModel(p, H_0, N_d) - B(:))./B(:);

maxIter = 200; lambda = 1e-2; tol = 1e-10;
p = params(:)';
r = resid(p);
for k = 1:maxIter
    G = jacobianFD(resid, p, free);
    dp = zeros(1,2);
    dp(free) = -((G'*G + lambda*diag(diag(G'*G)))\(G'*r))';
    p_new = p + dp;
    if any(p_new <= 0), lambda = lambda*10; continue; end
    r_new = resid(p_new);
    if norm(r_new) < norm(r)
        p = p_new; r = r_new; lambda = max(lambda/10, 1e-12);
        if norm(dp./p) < tol, break; end
    else
        lambda = lambda*10;
        if lambda > 1e12, break; end
    end
end
params = p;

G = jacobianFD(resid, p, free);
n = numel(B); np = sum(free);
s2 = sum(r.^2)/(n - np);
C = zeros(2);
C(free,free) = s2*inv(G'*G); %#ok<MINV>
info.cov = C;
info.sd = sqrt(diag(C))';
info.corr = NaN;
if ~fixMs, info.corr = C(1,2)/sqrt(C(1,1)*C(2,2)); end
info.rms = sqrt(mean(r.^2));
info.resid = r;
end

function G = jacobianFD(resid, p, free)
r0 = resid(p);
idx = find(free);
G = zeros(numel(r0), numel(idx));
for c = 1:numel(idx)
    j = idx(c);
    h = 1e-6*p(j);
    pp = p; pp(j) = p(j) + h;
    G(:,c) = (resid(pp) - r0)/h;
end
end
