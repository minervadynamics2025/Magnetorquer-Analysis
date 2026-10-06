function [ p_best ] = optimizeSimulatedAnnealingLangevinModel( N_d, B, H_0, M_s_fixed )
%OPTIMIZESIMULATEDANNEALINGLANGEVINMODEL Global search for the Langevin parameters.
% B         : volume-averaged flux density converted from the measurements [T]
% H_0       : applied magnetizing field H_app [A/m]
% N_d       : demagnetization factor
% M_s_fixed : (optional) fix M_s to this value [A/m] and search only for a
%
% Returns p_best = [M_s, a].
%
% Changes from the previous version
%  * mu0 = 4*pi*1e-7 everywhere (was 1.257e-6 here and 4*pi*1e-7 in the LM fit)
%  * the model is B = mu0*(H_eff + M) (the previous cost used mu0*(H_0 + M))
%  * relative residuals, so all measurement points have equal weight
%  * the search is done in log-parameter space (parameters span decades)

if nargin < 4, M_s_fixed = []; end
fixMs = ~isempty(M_s_fixed);

lb = log([1e4, 1e-1]);          % bounds of [M_s, a] (log space)
ub = log([4e6, 5e4]);
if fixMs, lb(1) = log(M_s_fixed); ub(1) = lb(1); end

cost = @(q) sum(((B - langevinRodModel(exp(q), H_0, N_d))./B).^2);

maxIter = 5000;
T = 1.0; T_min = 1e-6; alpha = 0.95;
q_current = lb + rand(1,2).*(ub - lb);
c_current = cost(q_current);
q_best = q_current; c_best = c_current;
for k = 1:maxIter
    step = 0.1*(ub - lb);
    q_new = min(max(q_current + step.*randn(1,2), lb), ub);
    c_new = cost(q_new);
    if c_new < c_current || rand < exp(-(c_new - c_current)/T)
        q_current = q_new; c_current = c_new;
    end
    if c_current < c_best
        q_best = q_current; c_best = c_current;
    end
    T = T*alpha;
    if T < T_min, break; end
end
p_best = exp(q_best);
end
