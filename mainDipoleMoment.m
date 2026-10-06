%% Relative permeability of a HyperCo 50A magnetorquer rod from teslameter data
% Revised analysis (manuscript revision R1):
%  1) core-loaded axisymmetric field solution -> transfer factor kappa and N_d
%     (replaces the air-core Green's-function factor r_s0*r_volume = 25.54)
%  2) Langevin fit: M_s free (effective parameters) and M_s fixed to datasheet
%  3) constant-mu_r reference fit, mu_r(H_eff) with propagated uncertainty
%  4) figures (PDF) and results_macros.tex for the manuscript
clear; clc; close all;

%% Geometry  >>> CONFIRM THE PLACEHOLDER VALUES <<<
geo.N      = 3210;      % number of turns
geo.L_wind = 0.150;     % winding length [m]
geo.L_rod  = 0.170;     % TOTAL rod length [m]                <-- PLACEHOLDER
geo.r0     = 0.005;     % rod radius [m]
geo.t_coil = 0.0026;    % radial build of the winding [m]     <-- PLACEHOLDER
geo.s      = 0.017;     % sensor distance from END OF WINDING [m]
d_copper   = 0.32e-3;   % wire diameter (AWG 28) [m]
B_sat      = 2.40;      % HyperCo 50A saturation induction, 24 kG (Carpenter) [T]
mu_ref     = 3000;      % starting mu_r for kappa; iterated to the fitted value
rho_cu     = 1.72e-8;   % copper resistivity at 20 C [ohm m]
R_coil     = 34.6;      % measured winding resistance [ohm]
probeCal   = 0.01;      % relative teslameter calibration uncertainty   <-- PLACEHOLDER
probeLen   = 1.0e-3;    % axial extent of the Hall element [m]            <-- PLACEHOLDER

mu0 = 4e-7*pi;

% Output folder for the manuscript (figures + numbers). Falls back to ./figs
manuscriptDir = ['C:\Users\anila\Dropbox\Minerva Dynamics Documents\', ...
    'Dizayn Dokümanları\Magnetorquer Design & Analysis\Relative Permeability Modeling'];
if ~exist(manuscriptDir, 'dir'), manuscriptDir = pwd; end
figDir = fullfile(manuscriptDir, 'figs');
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% Measurements
data = readmatrix('Tor-Dyn-25-10 Test Results.xlsx');
V      = data(:,1);              % applied voltage [V]
I      = data(:,2);              % measured current [A]
B_meas = data(:,3)/1000;         % teslameter reading [T]
H_0    = geo.N*I/geo.L_wind;     % applied magnetizing field H_app [A/m]

e_rod = (geo.L_rod - geo.L_wind)/2;
fprintf('Rod overhang %.1f mm, sensor-to-rod-tip %.1f mm\n', 1e3*e_rod, 1e3*(geo.s - e_rod));

%% 1) Core-loaded field solution: transfer factor and demagnetizing factor
H_0tmp = geo.N*I/geo.L_wind;
for it = 1:5                                      % evaluate kappa, N_d at the identified mu_r
    ref = getFieldMetrics(geo, mu_ref, I(end));
    mu_new = constantPermeabilityFit(ref.kappa*B_meas, H_0tmp, ref.Nd);
    if abs(mu_new/mu_ref - 1) < 2e-3, break; end
    mu_ref = mu_new;
end
fprintf('kappa and N_d evaluated at mu_r = %.0f\n', mu_ref);
air   = getFieldMetrics(geo, 1, I(end));          % coil-only (air-core) field
kappa = ref.kappa;
N_d   = ref.Nd;
Nd500 = getFieldMetrics(geo, 500, I(end)).Nd;
Nd1e5 = getFieldMetrics(geo, 1e5, I(end)).Nd;
kap_mu = arrayfun(@(m) getFieldMetrics(geo, m, I(end)).kappa, [1e3 1e4 1e5]);
kap_m1 = getFieldMetrics(geo, mu_ref, I(end), geo.s - 1e-3).kappa;
kap_p1 = getFieldMetrics(geo, mu_ref, I(end), geo.s + 1e-3).kappa;

kappaSpread  = 100*(max(kap_mu) - min(kap_mu))/kappa;      % % change, mu_r 1e3..1e5
dKappa       = 100*(kap_p1 - kap_m1)/(2*kappa);             % % per mm of standoff
kappaCoil    = air.kappa;                                   % air-core ratio
kappaM       = (ref.Bvol - air.Bvol)/(ref.Bs - air.Bs);     % magnetization-only ratio
MshareSensor = 100*(ref.Bs - air.Bs)/ref.Bs;
MshareInside = 100*(ref.Bvol - air.Bvol)/ref.Bvol;
fprintf('kappa = %.2f (coil-only %.2f, magnetization-only %.2f), N_d = %.5f\n', ...
    kappa, kappaCoil, kappaM, N_d);
fprintf('magnetization share: sensor %.1f %%, inside %.1f %%\n', MshareSensor, MshareInside);

% Measured field converted to the volume-averaged internal field
B = kappa*B_meas;
Bmax = mu0*H_0/N_d;                                         % mu_r -> infinity limit
if any(B >= Bmax)
    warning('Converted <B> exceeds the mu_r -> inf limit: check the geometry (L_rod, s).');
end

%% 2) Langevin fits (simulated annealing + Levenberg-Marquardt)
rng(1);
pA = optimizeSimulatedAnnealingLangevinModel(N_d, B, H_0);          % M_s free
[pA, infoA] = fitLMLangevin(B, H_0, N_d, pA);
M_s_sheet = B_sat/mu0;
pB = optimizeSimulatedAnnealingLangevinModel(N_d, B, H_0, M_s_sheet);% M_s fixed
[pB, infoB] = fitLMLangevin(B, H_0, N_d, pB, M_s_sheet);
[mu_c, rmsC] = constantPermeabilityFit(B, H_0, N_d);

% Coil geometry vs the measured winding resistance
L_wire = R_coil*(pi*d_copper^2/4)/rho_cu;
r_mean_R = L_wire/(2*pi*geo.N);                   % mean turn radius implied by R
r_mean_geo = geo.r0 + geo.t_coil/2;
geoAlt = geo; geoAlt.r_in = max(geo.r0, r_mean_R - geo.t_coil/2);
alt = getFieldMetrics(geoAlt, mu_ref, I(end));
mu_alt = constantPermeabilityFit(alt.kappa*B_meas, H_0, alt.Nd);
fprintf('winding mean radius: geometry %.1f mm, from R %.1f mm -> mu_r %.0f\n', 1e3*r_mean_geo, 1e3*r_mean_R, mu_alt);

fprintf('Free fit : M_s = %.3g +- %.2g A/m, a = %.3g +- %.2g A/m, corr = %.2f, RMS = %.2f %%\n', ...
    pA(1), infoA.sd(1), pA(2), infoA.sd(2), infoA.corr, 100*infoA.rms);
fprintf('  geometric ceiling H_app,max/N_d = %.3g A/m\n', max(H_0)/N_d);
fprintf('Fixed M_s: a = %.3g +- %.2g A/m, RMS = %.2f %%\n', pB(2), infoB.sd(2), 100*infoB.rms);
fprintf('Constant mu_r = %.0f, RMS = %.2f %%\n', mu_c, 100*rmsC);

%% 3) Adopted model (M_s fixed): M, H_eff, mu_r with uncertainty
[B_calc, M_f, H_f] = langevinRodModel(pB, H_0, N_d);
mu_r = B./(mu0*H_f);
gain = 1 + N_d*(mu_r - 1);                   % d ln(mu_r) / d ln(<B>)
sig  = max(infoB.rms, 0.01);
mu_lo = mu_r./(1 + gain*sig);
mu_hi = mu_r./max(1 - gain*sig, 1e-9);
m_d  = pi*geo.r0^2*geo.L_rod*M_f;            % core dipole moment [A m^2]

%% 4) Sensitivity to the sensor position (one field solution suffices)
d_tip = linspace(0.003, 0.012, 46);                      % sensor-to-rod-tip distance [m]
s_tip = geo.s - e_rod;
kap_d = ref.Bvol./interp1(ref.z, ref.Bax, geo.L_rod/2 + d_tip);
mur_d = arrayfun(@(k) constMuOrInf(k*B_meas, H_0, N_d, Bmax), kap_d);
k_far  = ref.Bvol/interp1(ref.z, ref.Bax, geo.L_rod/2 + s_tip + 5e-4);
k_near = ref.Bvol/interp1(ref.z, ref.Bax, geo.L_rod/2 + s_tip - 5e-4);
murPos = sort([constMuOrInf(k_far*B_meas, H_0, N_d, Bmax), constMuOrInf(k_near*B_meas, H_0, N_d, Bmax)]);
fprintf('mu_r for +-0.5 mm sensor position: %.0f .. %.0f\n', murPos);
m_ar = geo.L_rod/(2*geo.r0);                             % prolate ellipsoid, same aspect ratio
NdEll = (m_ar/sqrt(m_ar^2-1)*log(m_ar + sqrt(m_ar^2-1)) - 1)/(m_ar^2 - 1);
posTol = 20/(mean(gain)*abs(dKappa));                    % mm for 20 % error in mu_r

% Uncertainty budget for the constant mu_r (one input perturbed at a time)
murFor = @(fac, k, nd) constMuOrInf(fac*k*B_meas, H_0, nd, mu0*H_0/nd);
sBfit  = rmsC/sqrt(numel(B_meas) - 1);
budFit = (murFor(1 + sBfit, kappa, N_d) - murFor(1 - sBfit, kappa, N_d))/2/mu_c;
budCal = (murFor(1 + probeCal, kappa, N_d) - murFor(1 - probeCal, kappa, N_d))/2/mu_c;
zs0 = geo.L_wind/2 + geo.s;
BsAvg = mean(interp1(ref.z, ref.Bax, linspace(zs0 - probeLen/2, zs0 + probeLen/2, 41)));
budArea = murFor(1, ref.Bvol/BsAvg, N_d)/mu_c - 1;
budCoil = mu_alt/mu_c - 1;
fprintf('budget: fit %.1f%%, cal %.1f%%, area %.1f%%, coil %.1f%%\n', 100*[budFit budCal budArea budCoil]);

%% 5) Figures (IEEE single column, 3.4 in wide; print-safe palette)
C1 = [42 120 214]/255; C2 = [235 104 52]/255;            % blue, orange
INK2 = [82 81 78]/255; BAND = [230 228 221]/255; BAND2 = [241 240 236]/255;

% Fig. 3 - centreline field, total vs coil-only
zf = ref.z(ref.z <= 0.13); n = numel(zf);
zz = 100*[-fliplr(zf), zf];
f1 = newFig(1.85); ax = gca; hold on
patch(100*[-1 1 1 -1]*geo.L_rod/2, [1e-5 1e-5 1 1], BAND2, 'EdgeColor', 'none');
patch(100*[-1 1 1 -1]*geo.L_wind/2, [1e-5 1e-5 1 1], BAND, 'EdgeColor', 'none');
semilogy(zz, [fliplr(ref.Bax(1:n)), ref.Bax(1:n)], 'Color', C1, 'LineWidth', 1.5);
semilogy(zz, [fliplr(air.Bax(1:n)), air.Bax(1:n)], '--', 'Color', C2, 'LineWidth', 1.5);
zs_cm = 100*(geo.L_wind/2 + geo.s);
plot([zs_cm zs_cm], [1e-5 1], ':', 'Color', INK2);
plot(zs_cm, ref.Bs, 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', C1, 'LineWidth', 1.2);
text(zs_cm + 0.3, 2*ref.Bs, 'sensor', 'FontSize', 7, 'Color', INK2);
text(0, 2.5e-5, 'winding', 'FontSize', 7, 'Color', INK2, 'HorizontalAlignment', 'center');
set(ax, 'YScale', 'log'); xlim([-13 13]); ylim([1e-5 1]);
xlabel('z (cm)'); ylabel('B_z(0,z) (T)');
legend({'', '', 'coil + core (total)', 'coil only (air core)'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
styleAxes(ax); saveFig(f1, fullfile(figDir, 'rev_FEMCenterline.pdf'));

% Fig. 4 - converted measurements vs models
Hl = linspace(1, 1.05*max(H_0), 200)';
Bl = langevinRodModel(pB, Hl, N_d);
f2 = newFig(2.0); ax = gca; hold on
plot(Hl, mu0*Hl/N_d, ':', 'Color', INK2, 'LineWidth', 0.9);

text(Hl(end), mu0*Hl(end)/N_d, ' \mu_r\rightarrow\infty', 'FontSize', 7, 'Color', INK2);
text(Hl(end), mu0*1000*Hl(end)/(1 + N_d*999), ' \mu_r=1000', 'FontSize', 7, 'Color', INK2);
text(Hl(end), mu0*300*Hl(end)/(1 + N_d*299), ' \mu_r=300', 'FontSize', 7, 'Color', INK2);
hL = plot(Hl, Bl, 'Color', C1, 'LineWidth', 1.5);
plot(Hl, mu0*1000*Hl./(1 + N_d*999), '--', 'Color', [11 11 11]/255, 'LineWidth', 0.9);   % on top
plot(Hl, mu0*300*Hl./(1 + N_d*299), ':', 'Color', [11 11 11]/255, 'LineWidth', 0.9);
hM = errorbar(H_0, B, infoB.rms*B, 'o'); styleMarkers(hM, C2);
xlim([0 1.22*max(Hl)]); ylim([0 inf]);
xlabel('H_{app} (A/m)'); ylabel('\langle B_z\rangle (T)');
legend([hL hM], {'Langevin model (M_s fixed)', 'measured, converted with \kappa'}, 'Location', 'northwest');
styleAxes(ax); saveFig(f2, fullfile(figDir, 'rev_MagneticField.pdf'));

% Fig. 5 - sensitivity to the sensor position
f3 = newFig(2.45);
t = tiledlayout(f3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
a1 = nexttile(t); hold on
patch(1e3*(s_tip + [-5e-4 5e-4 5e-4 -5e-4]), [0 0 1 1]*1.1*max(kap_d), BAND, 'EdgeColor', 'none');
plot(1e3*d_tip, kap_d, 'Color', C1, 'LineWidth', 1.5);
ylabel('\kappa'); ylim([0.9*min(kap_d) 1.05*max(kap_d)]); set(a1, 'XTickLabel', []); styleAxes(a1);
a2 = nexttile(t); hold on
fin = isfinite(mur_d);
if any(~fin)
    patch(1e3*[min(d_tip(~fin)) max(d_tip(~fin)) max(d_tip(~fin)) min(d_tip(~fin))], [1 1 1e5 1e5], [246 213 199]/255, 'EdgeColor', 'none');
    text(1e3*mean(d_tip(~fin)), sqrt(min(mur_d(fin))*max(mur_d(fin))), {'data above', '\mu_r\rightarrow\infty limit'}, ...
        'HorizontalAlignment', 'center', 'FontSize', 6.5, 'Color', INK2);
end
patch(1e3*(s_tip + [-5e-4 5e-4 5e-4 -5e-4]), [1 1 1e5 1e5], BAND, 'EdgeColor', 'none');
plot(1e3*d_tip(fin), mur_d(fin), 'Color', C2, 'LineWidth', 1.5);
set(a2, 'YScale', 'log'); ylim([0.8*min(mur_d(fin)) 1.3*max(mur_d(fin))]);
xlabel('sensor distance from rod tip (mm)'); ylabel('fitted \mu_r'); styleAxes(a2);
saveFig(f3, fullfile(figDir, 'rev_Sensitivity.pdf'));

% Fig. 6 - relative permeability with propagated uncertainty
f4 = newFig(1.8); ax = gca; hold on
plot([0 1.05*max(H_f)], [mu_c mu_c], 'Color', C1, 'LineWidth', 1.2);
text(0.03*max(H_f), mu_c*1.05, sprintf('constant fit, \\mu_r = %.0f', mu_c), 'FontSize', 7, 'Color', INK2);
hE = errorbar(H_f, mu_r, mu_r - mu_lo, min(mu_hi - mu_r, 2*mu_r), 'o'); styleMarkers(hE, C2);
xlim([0 1.05*max(H_f)]); xlabel('H_{eff} (A/m)'); ylabel('\mu_r');
styleAxes(ax); saveFig(f4, fullfile(figDir, 'rev_RelativePermeability.pdf'));

%% 6) Numbers for the manuscript (results_macros.tex)
fmt = @(x, d) strrep(sprintf(['%.', num2str(d), 'f'], x), ',', '');
mac = {
 'kappaVal', fmt(kappa,1); 'NdVal', fmt(N_d,4); 'kappaSpread', fmt(kappaSpread,1);
 'dKappa', fmt(abs(dKappa),0); 'BpeakVal', fmt(ref.Bpeak,2); 'BvolVal', fmt(ref.Bvol,2);
 'BendVal', fmt(ref.Bend,3); 'BsVal', fmt(1e3*ref.Bs,2); 'BvolRatio', fmt(ref.Bvol/ref.Bend,0);
 'LrodVal', fmt(100*geo.L_rod,1); 'sTipVal', fmt(1e3*(geo.s - e_rod),0);
 'eRodVal', fmt(1e3*e_rod,0); 'tCoilVal', fmt(1e3*geo.t_coil,1);
 'MsFree', fmt(pA(1)/1e6,3); 'MsFreeSd', fmt(infoA.sd(1)/1e6,3); 'aFree', fmt(pA(2),1);
 'aFreeSd', fmt(infoA.sd(2),1); 'corrFree', fmt(infoA.corr,2); 'rmsFree', fmt(100*infoA.rms,1);
 'muMsFree', fmt(mu0*pA(1),2); 'MsSheet', fmt(M_s_sheet/1e6,2); 'aFix', fmt(pB(2),0);
 'aFixSd', fmt(infoB.sd(2),0); 'rmsFix', fmt(100*infoB.rms,1); 'murConst', fmt(mu_c,0);
 'rmsConst', fmt(100*rmsC,1); 'murMin', fmt(min(mu_r),0); 'murMax', fmt(max(mu_r),0);
 'gainMin', fmt(min(gain),0); 'gainMax', fmt(max(gain),0); 'HeffMax', fmt(max(H_f),0);
 'MMax', fmt(max(M_f)/1e6,3); 'MFrac', fmt(100*max(M_f)/M_s_sheet,0);
 'BmaxLim', fmt(max(Bmax),3); 'budFit', fmt(100*budFit,0); 'budCal', fmt(100*budCal,0);
 'budArea', fmt(100*abs(budArea),1); 'budCoil', fmt(100*abs(budCoil),1); 'probeCal', fmt(100*probeCal,0);
 'probeLen', fmt(1e3*probeLen,1); 'sBfit', fmt(100*sBfit,1); 'murRef', fmt(mu_ref,0); 'NdVar', fmt(100*(Nd500 - Nd1e5)/Nd500,1);
 'MprofMid', fmt(100*ref.Mmid,0); 'MprofEnd', fmt(100*ref.Mend,0);
 'rMeanGeo', fmt(1e3*r_mean_geo,1); 'rMeanR', fmt(1e3*r_mean_R,1); 'murCoilAlt', fmt(mu_alt,0);
 'NdCoilAlt', fmt(100*(alt.Nd/N_d - 1),0); 'kappaCoilAlt', fmt(100*abs(alt.kappa/kappa - 1),1);
 'NdMsVal', fmt(N_d*M_s_sheet/1e4,1); 'satRatio', fmt(N_d*M_s_sheet/max(H_0),0); 'murPosLo', fmt(murPos(1),0); 'murPosHi', infOr(murPos(2));
 'posTol', fmt(posTol,2); 'NdEll', fmt(NdEll,4); 'gammaVal', fmt(m_ar,0); 'mdipMax', fmt(max(m_d),2);
 'HappMax', thousands(max(H_0));
 'kappaCoil', fmt(kappaCoil,1); 'kappaM', fmt(kappaM,1); 'MshareSensor', fmt(MshareSensor,0);
 'MshareInside', fmt(MshareInside,0); 'kappaCoilErr', fmt(100*(kappa - kappaCoil)/kappa,0);
 'MceilFree', fmt(100*pA(1)/(max(H_0)/N_d),0); 'Mceil', fmt(max(H_0)/N_d/1e6,3)};
fid = fopen(fullfile(manuscriptDir, 'results_macros.tex'), 'w');
fprintf(fid, '%% Auto-generated by mainDipoleMoment.m -- do not edit by hand\n');
for i = 1:size(mac,1)
    fprintf(fid, '\\newcommand{\\%s}{%s}\n', mac{i,1}, mac{i,2});
end
fclose(fid);
fprintf('Wrote figures to %s and results_macros.tex\n', figDir);
