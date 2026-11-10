function out = ValidationXu2021(Npart,Ttau,fsbc,dtsub)
%VALIDATIONXU2021  Coupled body/Rayleigh multiple-scattering benchmark.
%
%   out = ValidationXu2021(Npart,Ttau,fsbc,dtsub) reproduces Figures 2-3 and
%   Table 1 of Xu, Margerin & Mikesell (2021), GJI 228, 1213-1236 (doi:10.1093/gji/ggab370): a point
%   source at 1 km depth in a homogeneous halfspace with point density
%   scatterers, run until the energy ratios between Rayleigh, P, Sv and Sh
%   waves settle on their equipartition values.
%
%     Npart  number of particles          (default 2e5)
%     Ttau   run length in units of tau_R (default 200)
%     fsbc   'modeconv' or 'mirror'       (default 'modeconv')
%     dtsub  propagation sub-step, s      (default 0.1)
%
%   Everything is in km and s. See also RAYLEIGHMODEL, TESTRAYLEIGHMODEL.

if nargin < 1 || isempty(Npart), Npart = 2e5;        end
if nargin < 2 || isempty(Ttau),  Ttau  = 200;        end
if nargin < 3 || isempty(fsbc),  fsbc  = 'modeconv'; end
if nargin < 4 || isempty(dtsub), dtsub = 0.1;        end

%% 1. medium and scatterers (Xu et al., Section 4)
alpha0 = 5.0;       % km/s
beta0  = 3.0;       % km/s
f0     = 1.0;       % Hz
V      = 0.5;       % km^3   (the paper prints m^3, see RayleighModel)
drho   = 0.2;       % delta_rho / rho0
n      = 1.0;       % km^-3
zmax   = 10.0;      % km, thickness of the scattering layer
z0     = 1.0;       % km, source depth

M = RayleighModel(alpha0,beta0,f0,V,drho,n,zmax);

fprintf('=== medium ===\n');
fprintf('  vP %.2f  vS %.2f  cR %.4f km/s   lambda_R %.3f km\n', ...
    alpha0,beta0,M.cR,M.cR/f0);
fprintf('=== mean free paths ===\n');
fprintf('  l_R  %7.2f km  (Xu et al. 35.8)   tau_R %6.2f s\n',M.lR,M.tauR);
fprintf('  l_P  %7.2f km   l_Sv %7.2f km   l_Sh %7.2f km\n',M.lB);
fprintf('=== theoretical equipartition (Eq. 17-19) ===\n');
fprintf('  E_R/E_P(1km)  %7.3f  (Xu et al. 8.310)\n',M.ratioRB(1));
fprintf('  E_R/E_Sv(1km) %7.3f  (Xu et al. 1.795)\n',M.ratioRB(2));
fprintf('  E_R/E_Sh(1km) %7.3f  (Xu et al. 1.795)\n',M.ratioRB(3));
fprintf('  E_Sv/E_P      %7.3f  = (vP/vS)^3\n',M.ratioSvP);

%% 2. Monte Carlo run
dtrec = 1.0;                        % s, recording interval (Xu et al. dt)
nsub  = max(1,round(dtrec/dtsub));
dtsub = dtrec/nsub;
T     = Ttau*M.tauR;
Nt    = floor(T/dtrec);
nbin  = round(zmax);                % 1 km depth bins

% source: impulsive P wave at 1 km depth, isotropic radiation
rng(1);
P.x    = [zeros(Npart,2) z0*ones(Npart,1)];
c  = 2*rand(Npart,1)-1;  ph = 2*pi*rand(Npart,1);  s = sqrt(1-c.^2);
P.dir  = [s.*cos(ph) s.*sin(ph) c];
P.mode = ones(Npart,1);

fprintf('=== Monte Carlo: %g particles, %g s = %.0f tau_R, dt = %g s, %s FSBC ===\n', ...
    Npart,Nt*dtrec,Nt*dtrec/M.tauR,dtsub,fsbc);

ER = zeros(Nt,1);                   % Rayleigh energy over the whole medium
EZ = zeros(Nt,nbin,3);              % body-wave energy density per 1 km bin
tic
for it = 1:Nt
    for is = 1:nsub
        P = propagateParticleRayleigh(M,P,dtsub,fsbc);
    end
    isR = P.mode == 4;
    ER(it) = sum(isR)/Npart;
    ib = floor(P.x(:,3))+1;
    ok = ~isR & ib >= 1 & ib <= nbin;
    EZ(it,:,:) = reshape(accumarray([ib(ok) P.mode(ok)],1,[nbin 3]),1,nbin,3)/Npart;
    if mod(it,max(1,round(Nt/10))) == 0
        fprintf('  t = %6.0f s (%3.0f%%)  E_R = %.4f  elapsed %.0f s\n', ...
            it*dtrec,100*it/Nt,ER(it),toc);
    end
end

%% 3. equipartition ratios
t    = (1:Nt).'*dtrec;
EB1  = squeeze(EZ(:,1,:));                     % shallowest 1 km layer
ratio = ER./EB1;
pl   = t > 0.75*t(end);                        % plateau = last quarter of the run
rpl  = mean(ER(pl))./mean(EB1(pl,:));          % ratio of the time-averaged energies
svp  = mean(sum(EZ(pl,:,2),2))/mean(sum(EZ(pl,:,1),2));   % over the whole layer
shp  = mean(sum(EZ(pl,:,3),2))/mean(sum(EZ(pl,:,1),2));

fprintf('=== simulated equipartition (mean over the last %.0f tau_R) ===\n', ...
    0.25*Nt*dtrec/M.tauR);
nm = {'E_R/E_P ','E_R/E_Sv','E_R/E_Sh'};
for i = 1:3
    fprintf('  %s %8.3f   theory %7.3f   misfit %5.1f%%\n', ...
        nm{i},rpl(i),M.ratioRB(i),100*abs(rpl(i)/M.ratioRB(i)-1));
end
fprintf('  E_Sv/E_P %8.3f   theory %7.3f   misfit %5.1f%%\n', ...
    svp,M.ratioSvP,100*abs(svp/M.ratioSvP-1));
fprintf('  E_Sh/E_P %8.3f   theory %7.3f   misfit %5.1f%%\n', ...
    shp,M.ratioSvP,100*abs(shp/M.ratioSvP-1));

out = struct('M',M,'t',t,'ER',ER,'EZ',EZ,'ratio',ratio,'plateau',rpl, ...
             'EsvEp',svp,'EshEp',shp,'Npart',Npart,'fsbc',fsbc);

%% 4. figures
zc = (0.5:1:nbin-0.5);
f1 = figure('Name','Energy ratios','Position',[100 100 700 500]);
semilogy(t/M.tauR,ratio,'LineWidth',1); hold on
set(gca,'ColorOrderIndex',1);
plot(t([1 end])/M.tauR,[M.ratioRB;M.ratioRB],'--','LineWidth',1.5);
xlabel('time (\tau_R)'); ylabel('energy ratio'); grid on
legend('E_R/E_P','E_R/E_{Sv}','E_R/E_{Sh}','theory','Location','best');
title(sprintf('Xu et al. (2021) Fig. 3 -- %g particles, %s FSBC',Npart,fsbc));

f2 = figure('Name','Depth profiles','Position',[820 100 900 320]);
snap = round([0.1 0.25 0.5 1]*Nt);
lbl  = arrayfun(@(k) sprintf('%.0f\\tau_R',t(k)/M.tauR),snap,'uni',0);
ttl  = {'P','Sv','Sh'};
for i = 1:3
    subplot(1,3,i);
    plot(zc,squeeze(EZ(snap,:,i)).','LineWidth',1);
    xlabel('depth (km)'); ylabel('energy density (km^{-1})');
    title(ttl{i}); grid on
    if i == 1, legend(lbl,'Location','best'); end
end

if ~usejava('desktop')
    exportgraphics(f1,'ValidationXu2021_ratios.png','Resolution',150);
    exportgraphics(f2,'ValidationXu2021_depth.png','Resolution',150);
    fprintf('figures written to ValidationXu2021_*.png\n');
end
end
