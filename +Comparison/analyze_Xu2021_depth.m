function out = analyze_Xu2021_depth(Npart, Ttau, dtsub)
%ANALYZE_XU2021_DEPTH  Is the E_Sv/E_P deficit a free-surface effect or a model error?
%
%   Comparison.test_Xu2021 reports E_Sv/E_P summed over the shallowest 10 km
%   and finds it ~8% below the theoretical (vP/vS)^3 = 4.630. That comparison
%   is ambiguous: (vP/vS)^3 is an equipartition relation of the DIFFUSIVE
%   INTERIOR, and with l_P = 147 km the whole 10 km layer sits inside one
%   tenth of a mean free path of the free surface, where P<->Sv conversion at
%   the surface reshapes the ratio. A uniform deficit would point at the
%   scattering model; a deficit confined to the shallow bins points at the
%   surface, and is physical.
%
%   This resolves it by binning FAR deeper than 10 km and reading the ratio
%   bin by bin, with a statistical error bar per bin. Same medium, same
%   source and same propagator as Comparison.test_Xu2021.
%
%   out.zc, out.ratioSvP, out.ratioShP, out.sigmaSvP, out.nP, ...
%
%   See also COMPARISON.TEST_XU2021, PROPAGATEPARTICLERAYLEIGH.

if nargin < 1 || isempty(Npart), Npart = 1e5; end
if nargin < 2 || isempty(Ttau),  Ttau  = 80;  end
if nargin < 3 || isempty(dtsub), dtsub = 0.1; end

alpha0 = 5.0; beta0 = 3.0; f0 = 1.0;
V = 0.5; drho = 0.2; n = 1.0; zmax = 10.0; z0 = 1.0;
M = RayleighModel(alpha0, beta0, f0, V, drho, n, zmax);

% Depth bins: fine where the surface matters, coarse where the interior is.
zEdges = [0:1:10, 12:2:30, 35:5:60, 70:10:120, 150, 200, 300];
nbin   = numel(zEdges)-1;
zc     = 0.5*(zEdges(1:end-1)+zEdges(2:end));

dtrec = 1.0; nsub = max(1,round(dtrec/dtsub)); dtsub = dtrec/nsub;
T = Ttau*M.tauR; Nt = floor(T/dtrec);

rng(1);
P.x = [zeros(Npart,2), z0*ones(Npart,1)];
c = 2*rand(Npart,1)-1; ph = 2*pi*rand(Npart,1); s = sqrt(1-c.^2);
P.dir = [s.*cos(ph), s.*sin(ph), c];
P.mode = ones(Npart,1);

fprintf('\n=== analyze_Xu2021_depth: %g particles, %.0f tau_R, %d depth bins to %g km ===\n', ...
    Npart, Ttau, nbin, zEdges(end));
fprintf('l_P = %.1f km, l_Sv = %.1f km, l_R = %.1f km\n', M.lB(1), M.lB(2), M.lR);

EZ = zeros(Nt, nbin, 3);
for it = 1:Nt
    for is = 1:nsub
        P = propagateParticleRayleigh(M, P, dtsub, 'modeconv');
    end
    isB = P.mode < 4;
    ib = discretize(P.x(isB,3), zEdges);
    md = P.mode(isB);
    ok = ~isnan(ib);
    EZ(it,:,:) = reshape(accumarray([ib(ok), md(ok)], 1, [nbin 3]), 1, nbin, 3);
    if mod(it, max(1,round(Nt/10))) == 0
        fprintf('  %3.0f%%\n', 100*it/Nt);
    end
end

%% Plateau averages, per depth bin
t  = (1:Nt).'*dtrec;
pl = t > 0.70*t(end);
nIndep = max(1, 0.30*Ttau);

nP  = squeeze(mean(EZ(pl,:,1),1)).';
nSv = squeeze(mean(EZ(pl,:,2),1)).';
nSh = squeeze(mean(EZ(pl,:,3),1)).';

ratioSvP = nSv ./ max(nP,eps);
ratioShP = nSh ./ max(nP,eps);
relErr   = @(a,b) sqrt((1./max(a,1) + 1./max(b,1))/nIndep);
sigmaSvP = ratioSvP .* relErr(nSv,nP);
sigmaShP = ratioShP .* relErr(nSh,nP);

th = M.ratioSvP;
fprintf('\n  %-14s %8s %10s %10s %10s %8s\n', ...
    'depth (km)','N_P','E_Sv/E_P','+/- 1sig','dev vs 4.630','n-sigma');
for k = 1:nbin
    if nP(k) < 5, continue; end
    dev = ratioSvP(k)/th - 1;
    fprintf('  %6.1f-%-6.1f %8.0f %10.3f %10.3f %9.1f%% %8.1f\n', ...
        zEdges(k), zEdges(k+1), nP(k), ratioSvP(k), sigmaSvP(k), 100*dev, ...
        abs(ratioSvP(k)-th)/max(sigmaSvP(k),eps));
end

%% Verdict
shallow = zc <= 10;
deep    = zc >= 60 & nP.' >= 20;
wMean = @(v,s,m) sum(v(m)./s(m).^2)/sum(1./s(m).^2);
rShallow = wMean(ratioSvP, sigmaSvP, shallow.');
rDeep    = wMean(ratioSvP, sigmaSvP, deep.');
sDeep    = 1/sqrt(sum(1./sigmaSvP(deep.').^2));
sShallow = 1/sqrt(sum(1./sigmaSvP(shallow.').^2));

fprintf('\n--- verdict ---\n');
fprintf('  theory (vP/vS)^3                      = %.4f\n', th);
fprintf('  weighted mean,  0-10 km (shallow)     = %.4f +/- %.4f  (%+.1f%%, %.1f sigma)\n', ...
    rShallow, sShallow, 100*(rShallow/th-1), abs(rShallow-th)/sShallow);
fprintf('  weighted mean, >=60 km (interior)     = %.4f +/- %.4f  (%+.1f%%, %.1f sigma)\n', ...
    rDeep, sDeep, 100*(rDeep/th-1), abs(rDeep-th)/sDeep);
if abs(rDeep-th)/sDeep < 3 && abs(rShallow-th)/sShallow > 3
    fprintf('  => deficit is CONFINED TO THE SURFACE: physical, not a model error.\n');
elseif abs(rDeep-th)/sDeep > 3
    fprintf('  => deficit PERSISTS in the interior: points at the scattering model.\n');
else
    fprintf('  => inconclusive at this statistics; raise Npart or Ttau.\n');
end

%% Figure
f = figure('Name','E_Sv/E_P vs depth','Color','w','Position',[100 100 860 520]);
keep = nP >= 5;
errorbar(zc(keep), ratioSvP(keep), sigmaSvP(keep), 'o-','LineWidth',1.2,'MarkerSize',4); hold on
errorbar(zc(keep), ratioShP(keep), sigmaShP(keep), 's-','LineWidth',1.2,'MarkerSize',4);
yline(th,'k--','(v_P/v_S)^3 = 4.630','LineWidth',1.5);
xline(10,'r:','10 km (test\_Xu2021 layer)','LineWidth',1.2);
xline(M.lB(1),'m:','l_P','LineWidth',1.2);
set(gca,'XScale','log'); grid on
xlabel('depth (km)'); ylabel('energy ratio');
legend('E_{Sv}/E_P','E_{Sh}/E_P','theory','Location','best');
title(sprintf('Equipartition vs depth -- %g particles, %.0f \\tau_R', Npart, Ttau));
subtitle('the 10 km layer test\_Xu2021 sums over is well inside one l_P of the surface');
%if ~usejava('desktop')
    exportgraphics(f,'analyze_Xu2021_depth.png','Resolution',150);
    fprintf('figure written to analyze_Xu2021_depth.png\n');
%end

out = struct('M',M,'zEdges',zEdges,'zc',zc,'nP',nP,'nSv',nSv,'nSh',nSh, ...
             'ratioSvP',ratioSvP,'ratioShP',ratioShP, ...
             'sigmaSvP',sigmaSvP,'sigmaShP',sigmaShP, ...
             'rShallow',rShallow,'rDeep',rDeep,'theory',th);
end
