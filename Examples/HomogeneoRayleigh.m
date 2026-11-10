function R = HomogeneoRayleigh(bodyFile)
%HOMOGENEORAYLEIGH  Rayleigh channel for the cylinder of Homogeneo.m.
%
%   R = HomogeneoRayleigh() runs the Rayleigh-wave channel on the same
%   geometry, source, material and observation grid as Homogeneo.m and
%   plots the result.
%
%   R = HomogeneoRayleigh(bodyFile) also loads the body-wave field from a
%   previous radiativeTransfer run (e.g. 'HomogeneoUniform.mat') and
%   overlays the two, which is the point of the exercise: the Rayleigh
%   front must arrive at the rim later than S and must stay within about
%   two Rayleigh wavelengths of the free surface.
%
%   See also RAYLEIGHHOMOGENEOUS, TESTRAYLEIGHHOMOGENEOUS.

if nargin < 1, bodyFile = ''; end
base = '/Users/lac/Nextcloud/ArtigoManda/RTE/Elastico';

%% setup copied from Homogeneo.m
geometry = struct('dimension',3,'frame','cylindrical');
geometry.bnd(1) = struct('dir',3,'val',0,'axis',0);
geometry.bnd(2) = struct('dir',3,'val',-0.23,'axis',0);
geometry.bnd(3) = struct('dir',4,'val',0.055,'axis',3);

material = struct('vp',4814.4,'vs',2962.1,'rho',2210.1,'Frequency',500e3);

ic = load(fullfile(base,'ICrteManda2.mat'));
source = struct('numberParticles',5e6,'type','point','position',[1e-9 0 0], ...
                'direction','uniform','radial',griddedInterpolant(ic.r,ic.amp), ...
                'polarization','P');

xx = [0 2.5e-3:5e-4:5.5e-2];
observation = struct('x',xx,'y',[-pi pi],'z',-0.23:5e-4:0, ...
                     'directions',[0 pi],'time',0:0.5e-6:2e-4);

%% Rayleigh channel
tic
R = rayleighHomogeneous(geometry,source,material,observation, ...
        struct('numberParticles',4e5));
fprintf('Rayleigh channel: %.0f s\n',toc);

fprintf('cR = %.1f m/s, lambda_R = %.3f mm\n',R.cR,1e3*R.lambdaR);
fprintf('source partition  P %.1f%%  S %.1f%%  R %.1f%%\n',100*R.partition);
fprintf('edge partition    round %.0f%%  back %.0f%%  radiated %.0f%%\n',100*R.edge);
fprintf('rim arrival: R %.2f us, S %.2f us, P %.2f us\n', ...
    1e6*R.rmax/R.cR,1e6*R.rmax/material.vs,1e6*R.rmax/material.vp);
fprintf('Rayleigh energy left at %.0f us: %.1f%% of the source (%.1f%% radiated)\n', ...
    1e6*R.t(end),100*R.energy(end),100*R.radiated(end));

%% figures
skin = R.z > -2*R.lambdaR;                % the first two Rayleigh wavelengths

f1 = figure('Position',[80 80 1000 380]);
subplot(1,2,1)
plot(1e6*R.t,R.energy,'LineWidth',1.2); hold on
plot(1e6*R.t,R.radiated,'LineWidth',1.2);
plot(1e6*R.t,R.energy+R.radiated,'k--');
xlabel('t (\mus)'); ylabel('energy / source energy'); grid on
legend('Rayleigh','radiated at the edges','total','Location','east');
title('Rayleigh energy budget');

subplot(1,2,2)
imagesc(1e3*R.x,1e6*R.t,squeeze(sum(R.energyDensity(:,skin,:),2)).');
set(gca,'YDir','normal'); hold on
h1 = plot(1e3*R.cR*R.t,1e6*R.t,'r--','LineWidth',1);
h2 = plot(1e3*material.vs*R.t,1e6*R.t,'w:','LineWidth',1);
xlabel('r (mm)'); ylabel('t (\mus)'); ylim([0 60]);
legend([h1 h2],{'c_R t','v_S t'},'Location','southeast');
title('surface energy, top 2\lambda_R'); colormap(gca,'hot');

f2 = figure('Position',[80 500 1000 380]);
its = [8 20 41 61];
for k = 1:4
    subplot(1,4,k)
    imagesc(1e3*R.x,1e3*R.z,squeeze(R.energyDensity(:,:,its(k))).');
    set(gca,'YDir','normal'); ylim(1e3*[-0.04 0]);
    xlabel('r (mm)'); if k==1, ylabel('z (mm)'); end
    title(sprintf('t = %.1f \\mus',1e6*R.t(its(k))));
    colormap(gca,'hot');
end
sgtitle('Rayleigh energy density: the sheet hugs z = 0 with the r_1^2+r_2^2 profile');

%% optional overlay with a body-wave run
if ~isempty(bodyFile)
    S = load(bodyFile,'obs');
    EB = S.obs.energyDensity;
    % The body-wave run normalises the WHOLE source energy into P and S,
    % because it knows nothing about the Rayleigh channel. Under the
    % Miller-Pursey partition only fB of the source goes into body waves,
    % so the body curves are rescaled to put both on the same source.
    fB = sum(R.partition(1:2));
    nt = min(numel(R.t),size(EB,3));
    eP = fB*squeeze(sum(EB(:,skin,1:nt,1),[1 2]));
    eS = fB*squeeze(sum(EB(:,skin,1:nt,2),[1 2]));
    eR = squeeze(sum(R.energyDensity(:,skin,1:nt),[1 2]));
    tt = 1e6*R.t(1:nt);

    f3 = figure('Position',[80 80 1050 400]);
    subplot(1,2,1)
    semilogy(tt,eP,'LineWidth',1.2); hold on
    semilogy(tt,eS,'LineWidth',1.2);
    semilogy(tt,eR,'LineWidth',1.6);
    xline(1e6*R.rmax/material.vp,':','P na borda');
    xline(1e6*R.rmax/material.vs,':','S na borda');
    xline(1e6*R.rmax/R.cR,'r--','R na borda');
    xlabel('t (\mus)'); ylabel('energy density, top 2\lambda_R'); grid on
    xlim([0 60]); legend('P','S','Rayleigh','Location','southwest');
    title(sprintf('surface skin (body rescaled by f_B = %.3f)',fB));

    % where each mode sits in depth, averaged over the first rim transit
    subplot(1,2,2)
    w = R.t(1:nt) <= R.rmax/R.cR;
    semilogx(fB*squeeze(sum(EB(:,:,w,1),[1 3])),1e3*R.z,'LineWidth',1.2); hold on
    semilogx(fB*squeeze(sum(EB(:,:,w,2),[1 3])),1e3*R.z,'LineWidth',1.2);
    semilogx(squeeze(sum(R.energyDensity(:,:,w),[1 3])),1e3*R.z,'LineWidth',1.6);
    yline(-1e3*R.lambdaR,'k:','\lambda_R');
    ylim(1e3*[-0.05 0]); xlabel('energy density'); ylabel('z (mm)'); grid on
    legend('P','S','Rayleigh','Location','southwest');
    title('depth reach before the first rim arrival');
end

if ~usejava('desktop')
    exportgraphics(f1,'HomogeneoRayleigh_budget.png','Resolution',150);
    exportgraphics(f2,'HomogeneoRayleigh_snapshots.png','Resolution',150);
    if exist('f3','var')
        exportgraphics(f3,'HomogeneoRayleigh_overlay.png','Resolution',150);
    end
    fprintf('figures written to HomogeneoRayleigh_*.png\n');
end
end
