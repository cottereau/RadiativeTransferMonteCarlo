function [allPass, out] = test_Xu2021(Npart, Ttau, fsbc, dtsub)
%TEST_XU2021  Literature comparison benchmark against Xu et al. (2021).
%
%   [allPass, out] = Comparison.test_Xu2021(Npart, Ttau, fsbc, dtsub) reproduces
%   the equipartition benchmark of Xu, Margerin & Mikesell (2021), GJI 228,
%   1213-1236, doi:10.1093/gji/ggab370 (Figures 2-3 and Table 1): a point source at 1 km depth in a
%   homogeneous halfspace with point density scatterers, coupled to the
%   free surface via Rayleigh modes.
%
%   PASS/FAIL CRITERIA:
%     At late times, the energy ratios settle onto their theoretical
%     equipartition values -- but each ratio is read WHERE ITS THEORY HOLDS:
%
%       E_R / E_P, E_R / E_Sv, E_R / E_Sh   in the shallowest 1 km bin.
%         Eq. 17-19 of Xu et al. (2021) are derived for the surface layer,
%         so that is where they apply:
%           E_R / E_P  = pi * alpha0^3 / (omega * cR * UR)
%           E_R / E_Sv = pi * beta0^3  / (omega * cR * UR)
%
%       E_Sv / E_P, E_Sh / E_P              in the DIFFUSIVE INTERIOR,
%         z >= l_P. (alpha0/beta0)^3 is a bulk relation and does NOT hold
%         near a free surface, where P<->Sv conversion on reflection biases
%         the mix. Measured here (Comparison.analyze_Xu2021_depth, 3e5
%         particles), the deficit is -8.4% at 0-10 km and decays monotonically
%         to zero over one P mean free path:
%
%             0- 10 km   4.240  (-8.4%,  9.8 sigma from theory)
%            20- 30 km   4.39   (-5%)
%            50- 60 km   4.51   (-2.5%)
%           100-120 km   4.60   (-0.5%)
%           150-200 km   4.629  (-0.0%)
%           >= 60 km     4.6244 +/- 0.0096  (-0.1%, 0.5 sigma)  <-- theory
%
%         Summing over the shallowest 10 km, as this test used to, compares a
%         surface-biased average against an interior relation and fails for a
%         physical reason that has nothing to do with the code.
%
%   INPUTS (all optional):
%     Npart - Number of Monte Carlo particles (default 5e4 for fast test, 2e5 paper)
%     Ttau  - Simulation duration in units of tau_R (default 100)
%     fsbc  - Free surface boundary condition: 'modeconv' (default) or 'mirror'
%     dtsub - Time sub-step, s (default 0.1)
%
%   OUTPUTS:
%     allPass - True if all equipartition misfits are within statistical tolerance
%     out     - Detailed results struct (M, t, ER, EZ, ratio, plateau, misfit)
%
%   See also RAYLEIGHMODEL, PROPAGATEPARTICLERAYLEIGH, COMPARISON.TEST_SATO1994_V1.

if nargin < 1 || isempty(Npart), Npart = 1e5;        end
if nargin < 2 || isempty(Ttau),  Ttau  = 80;         end
if nargin < 3 || isempty(fsbc),  fsbc  = 'modeconv'; end
if nargin < 4 || isempty(dtsub), dtsub = 0.1;        end

fprintf('\n=== Comparison.test_Xu2021 (Xu, Margerin & Mikesell 2021) ===\n');

%% 1. Medium and scatterers (Xu et al. Section 4)
alpha0 = 5.0;       % km/s
beta0  = 3.0;       % km/s
f0     = 1.0;       % Hz
V      = 0.5;       % km^3
drho   = 0.2;       % delta_rho / rho0
n      = 1.0;       % km^-3
zmax   = 10.0;      % km
z0     = 1.0;       % km

M = RayleighModel(alpha0, beta0, f0, V, drho, n, zmax);

fprintf('Medium properties:\n');
fprintf('  vP = %.2f km/s, vS = %.2f km/s, cR = %.4f km/s, lambda_R = %.3f km\n', ...
    alpha0, beta0, M.cR, M.cR/f0);
fprintf('Mean free paths:\n');
fprintf('  l_R = %7.2f km (tau_R = %6.2f s), l_P = %7.2f km, l_S = %7.2f km\n', ...
    M.lR, M.tauR, M.lB(1), M.lB(2));
fprintf('Theoretical equipartition ratios (Eq. 17-19):\n');
fprintf('  E_R / E_P(1km)  = %7.3f\n', M.ratioRB(1));
fprintf('  E_R / E_Sv(1km) = %7.3f\n', M.ratioRB(2));
fprintf('  E_R / E_Sh(1km) = %7.3f\n', M.ratioRB(3));
fprintf('  E_Sv / E_P      = %7.3f = (vP/vS)^3\n\n', M.ratioSvP);

%% 2. Monte Carlo simulation
dtrec = 1.0;
nsub  = max(1, round(dtrec / dtsub));
dtsub = dtrec / nsub;
T     = Ttau * M.tauR;
Nt    = floor(T / dtrec);

% Depth bins: 1 km near the surface, where Eq. 17-19 is read, then coarsening
% out past l_P, where (alpha0/beta0)^3 is read. Scatterers fill the whole
% halfspace (see propagateParticleRayleigh), so the deep bins are populated.
zDeep  = M.lB(1);                                   % l_P: start of the interior
zEdges = unique([0:1:zmax, 15:5:60, 70:10:140, zDeep, 200, 300, 500]);
nbin   = numel(zEdges) - 1;
zc     = 0.5*(zEdges(1:end-1) + zEdges(2:end));
isDeep = zc >= zDeep;
fprintf('Depth binning: %d bins to %g km; interior region z >= l_P = %.1f km (%d bins)\n\n', ...
    nbin, zEdges(end), zDeep, nnz(isDeep));

rng(1);
P.x = [zeros(Npart, 2), z0 * ones(Npart, 1)];
c = 2 * rand(Npart, 1) - 1;
ph = 2 * pi * rand(Npart, 1);
s = sqrt(1 - c.^2);
P.dir  = [s.*cos(ph), s.*sin(ph), c];
P.mode = ones(Npart, 1);

fprintf('Running Monte Carlo: %g particles, %g s = %.0f tau_R (dt = %g s)...\n', ...
    Npart, Nt * dtrec, Nt * dtrec / M.tauR, dtsub);

ER = zeros(Nt, 1);
EZ = zeros(Nt, nbin, 3);

for it = 1:Nt
    for is = 1:nsub
        P = propagateParticleRayleigh(M, P, dtsub, fsbc);
    end
    if mod(it, max(1,round(Nt/10))) == 0
        fprintf('  %3.0f%%  (t = %6.0f s, E_R = %.4f)\n', 100*it/Nt, it*dtrec, sum(P.mode==4)/Npart);
    end
    isR = P.mode == 4;
    ER(it) = sum(isR) / Npart;
    ib = discretize(P.x(:, 3), zEdges);
    ok = ~isR & ~isnan(ib);
    EZ(it, :, :) = reshape(accumarray([ib(ok), P.mode(ok)], 1, [nbin, 3]), 1, nbin, 3) / Npart;
end

%% 3. Equipartition verification
t     = (1:Nt).' * dtrec;
EB1   = squeeze(EZ(:, 1, :));
ratio = ER ./ EB1;
pl    = t > 0.70 * t(end);

rpl   = mean(ER(pl)) ./ mean(EB1(pl, :));

% Bulk ratios from the diffusive interior only (see the header).
deepP  = mean(sum(EZ(pl, isDeep, 1), 2));
deepSv = mean(sum(EZ(pl, isDeep, 2), 2));
deepSh = mean(sum(EZ(pl, isDeep, 3), 2));
svp = deepSv / deepP;
shp = deepSh / deepP;

% The shallow value is kept as a diagnostic: it is expected to sit BELOW
% theory, and by how much is a physical statement about the free surface.
shallowP  = mean(sum(EZ(pl, zc <= zmax, 1), 2));
shallowSv = mean(sum(EZ(pl, zc <= zmax, 2), 2));
svp_shallow = shallowSv / shallowP;

misfit_ER_EP  = abs(rpl(1) / M.ratioRB(1) - 1);
misfit_ER_ESv = abs(rpl(2) / M.ratioRB(2) - 1);
misfit_ER_ESh = abs(rpl(3) / M.ratioRB(3) - 1);
misfit_ESv_EP = abs(svp / M.ratioSvP - 1);
misfit_ESh_EP = abs(shp / M.ratioSvP - 1);

fprintf('\nSimulated Equipartition Ratios (plateau average):\n');
fprintf('  E_R/E_P:   sim = %7.3f | theory = %7.3f | misfit = %5.1f%%\n', ...
    rpl(1), M.ratioRB(1), 100 * misfit_ER_EP);
fprintf('  E_R/E_Sv:  sim = %7.3f | theory = %7.3f | misfit = %5.1f%%\n', ...
    rpl(2), M.ratioRB(2), 100 * misfit_ER_ESv);
fprintf('  E_R/E_Sh:  sim = %7.3f | theory = %7.3f | misfit = %5.1f%%\n', ...
    rpl(3), M.ratioRB(3), 100 * misfit_ER_ESh);
fprintf('  E_Sv/E_P:  sim = %7.3f | theory = %7.3f | misfit = %5.1f%%   (z >= %.0f km)\n', ...
    svp, M.ratioSvP, 100 * misfit_ESv_EP, zDeep);
fprintf('  E_Sh/E_P:  sim = %7.3f | theory = %7.3f | misfit = %5.1f%%   (z >= %.0f km)\n', ...
    shp, M.ratioSvP, 100 * misfit_ESh_EP, zDeep);
fprintf('  [diagnostic] E_Sv/E_P over 0-%g km = %.3f (%+.1f%% -- free-surface\n', ...
    zmax, svp_shallow, 100*(svp_shallow/M.ratioSvP - 1));
fprintf('               bias, not an error; recovers over one l_P)\n');

%% 4. Statistical tolerance, per ratio, from the counts that actually entered
% the plateau average -- not a fixed number.
%
% Each ratio is a ratio of two particle counts, so its relative error is
%
%     sigma_r / r = sqrt( (1/N_A + 1/N_B) / n_indep )
%
% N_A and N_B are the typical INSTANTANEOUS counts in the plateau (Poisson,
% var = mean) and n_indep is how many independent snapshots the plateau
% average contains. Successive samples are strongly correlated -- they are
% the same particles one second later -- so the plateau's Nt samples do NOT
% all count. Particles decorrelate on the Rayleigh mean free time tau_R, and
% the plateau spans 0.30*T = 0.30*Ttau*tau_R, hence n_indep ~ 0.30*Ttau.
%
% So the tolerance scales as 1/sqrt(Npart * Ttau), which is the point: raising
% either tightens the test automatically instead of leaving a stale constant.
% Treating A and B as independent is conservative -- they are anticorrelated,
% a particle that turns Rayleigh leaves the body pool -- so the real error is
% somewhat smaller than this.
%
% WHY THIS MATTERS HERE: E_P in the shallowest 1 km bin is the scarcest
% channel in the whole run. At Npart = 1e5 / Ttau = 80 it holds only a few
% tens of particles, so a fixed 15% tolerance was testing almost nothing on
% that ratio while being needlessly tight on E_R/E_Sv.
kSigma   = 4;
tolFloor = 0.02;   % below this, 1 km binning and dtsub dominate over counting
nIndep   = max(1, 0.30 * Ttau);

NR    = Npart * mean(ER(pl));                       % Rayleigh particles
NB1   = Npart * mean(EB1(pl,:), 1);                 % P, Sv, Sh in the top 1 km
NBall = Npart * [deepP, deepSv, deepSh];            % P, Sv, Sh in the interior

relErr = @(a,b) sqrt((1./max(a,1) + 1./max(b,1)) / nIndep);
sigma_ratio = [relErr(NR,      NB1(1)), ...
               relErr(NR,      NB1(2)), ...
               relErr(NR,      NB1(3)), ...
               relErr(NBall(2), NBall(1)), ...
               relErr(NBall(3), NBall(1))];
tol_ratio = max(kSigma * sigma_ratio, tolFloor);

misfits = [misfit_ER_EP, misfit_ER_ESv, misfit_ER_ESh, misfit_ESv_EP, misfit_ESh_EP];
rname   = {'E_R/E_P ','E_R/E_Sv','E_R/E_Sh','E_Sv/E_P','E_Sh/E_P'};
counts  = [NB1(1), NB1(2), NB1(3), NBall(1), NBall(1)];
fprintf('  (E_R/* read in the 0-1 km bin; E_S*/E_P read at z >= %.0f km)\n', zDeep);

fprintf('\nStatistical tolerance (%d-sigma, n_indep = %.0f from Ttau = %g):\n', ...
    kSigma, nIndep, Ttau);
fprintf('  %-9s %10s %10s %10s %8s   %s\n','ratio','misfit','1-sigma','tol','N_denom','');
for k = 1:5
    if misfits(k) < tol_ratio(k), tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  %-9s %9.2f%% %9.2f%% %9.2f%% %8.0f   [%s]\n', rname{k}, ...
        100*misfits(k), 100*sigma_ratio(k), 100*tol_ratio(k), counts(k), tag);
end

allPass = all(misfits < tol_ratio);

% A tolerance this loose means the run cannot discriminate a correct model
% from a wrong one; say so rather than let a weak PASS look strong.
if max(tol_ratio) > 0.15
    fprintf(['  WARNING: tolerance reaches %.1f%% -- this run is under-powered.\n' ...
             '           A PASS does not discriminate. Raise Npart or Ttau\n' ...
             '           (tolerance falls as 1/sqrt(Npart*Ttau)).\n'], ...
             100*max(tol_ratio));
end

if allPass
    fprintf('>>> Comparison.test_Xu2021: ALL EQUIPARTITION RATIOS PASS! <<<\n\n');
else
    fprintf('>>> Comparison.test_Xu2021: SOME RATIOS DEVIATED <<<\n\n');
end

%% 5. Figures
f1 = figure('Name','Xu2021: energy ratios','Color','w','Position',[100 100 760 520]);
semilogy(t/M.tauR, ratio, 'LineWidth', 1); hold on
set(gca,'ColorOrderIndex',1);
plot(t([1 end])/M.tauR, [M.ratioRB; M.ratioRB], '--', 'LineWidth', 1.5);
% Tolerance band: what the PASS/FAIL above is actually asking of each ratio.
set(gca,'ColorOrderIndex',1);
co = get(gca,'ColorOrder');
for k = 1:3
    band = M.ratioRB(k) * [1-tol_ratio(k), 1+tol_ratio(k)];
    fill([t(1) t(end) t(end) t(1)]/M.tauR, [band(1) band(1) band(2) band(2)], ...
         co(k,:), 'FaceAlpha', 0.12, 'EdgeColor','none', 'HandleVisibility','off');
end
xline(0.70*t(end)/M.tauR, 'k:', 'plateau', 'LineWidth', 1.2, ...
      'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
xlabel('time (\tau_R)'); ylabel('energy ratio'); grid on
legend('E_R/E_P','E_R/E_{Sv}','E_R/E_{Sh}','theory','Location','best');
title(sprintf('Xu et al. (2021) Fig. 3 -- %g particles, %.0f \\tau_R, %s FSBC', ...
    Npart, Ttau, fsbc));
subtitle('shaded = statistical tolerance band actually tested');

f2 = figure('Name','Xu2021: depth profiles','Color','w','Position',[880 100 960 620]);
snap = unique(max(1, round([0.1 0.25 0.5 1]*Nt)));
lbl  = arrayfun(@(k) sprintf('%.0f\\tau_R', t(k)/M.tauR), snap, 'uni', 0);
ttl  = {'P','Sv','Sh'};
binW = diff(zEdges);                    % bins are not uniform: divide it out
for i = 1:3
    subplot(2,2,i);
    plot(zc, squeeze(EZ(snap,:,i)).'./binW.', 'LineWidth', 1);
    set(gca,'XScale','log');
    xlabel('depth (km)'); ylabel('energy density (km^{-1})');
    title(ttl{i}); grid on
    if i == 1, legend(lbl, 'Location','best'); end
end

% The panel that justifies where each ratio is read.
subplot(2,2,4);
nPz  = squeeze(mean(EZ(pl,:,1),1)).';
nSvz = squeeze(mean(EZ(pl,:,2),1)).';
nShz = squeeze(mean(EZ(pl,:,3),1)).';
good = Npart*nPz >= 5;
rz   = nSvz./max(nPz,eps);
sz   = rz .* sqrt((1./max(Npart*nSvz,1) + 1./max(Npart*nPz,1))/nIndep);
errorbar(zc(good), rz(good), sz(good), 'o-', 'LineWidth', 1.1, 'MarkerSize', 4); hold on
yline(M.ratioSvP, 'k--', '(v_P/v_S)^3', 'LineWidth', 1.5);
xline(zmax,  'r:', 'surface layer', 'LineWidth', 1.2);
xline(zDeep, 'm:', 'l_P', 'LineWidth', 1.2);
set(gca,'XScale','log'); grid on
xlabel('depth (km)'); ylabel('E_{Sv}/E_P');
title('where the bulk ratio is read');
sgtitle('Xu et al. (2021) Fig. 2 -- depth profiles, and the equipartition depth scale');

if ~usejava('desktop')
    exportgraphics(f1,'test_Xu2021_ratios.png','Resolution',150);
    exportgraphics(f2,'test_Xu2021_depth.png','Resolution',150);
    fprintf('figures written to test_Xu2021_*.png\n');
end

out = struct('M', M, 't', t, 'ER', ER, 'EZ', EZ, 'ratio', ratio, ...
             'plateau', rpl, 'EsvEp', svp, 'EshEp', shp, ...
             'misfits', misfits, 'sigma_ratio', sigma_ratio, ...
             'tol_ratio', tol_ratio, 'nIndep', nIndep, ...
             'counts', struct('NR',NR,'NB1',NB1,'NBdeep',NBall), ...
             'zEdges', zEdges, 'zc', zc, 'zDeep', zDeep, 'isDeep', isDeep, ...
             'svp_shallow', svp_shallow, ...
             'allPass', allPass);

end
