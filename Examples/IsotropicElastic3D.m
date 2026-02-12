close all
clear all
clc

titlecase = '3D elastic case with isotropic scattering';
disp(['Testing ' titlecase ' ...']);

geometry = struct( 'dimension', 3 );

source = struct( 'numberParticles', 1e6, ...
    'position', [0 0 0], ...
    'polarization', 'P', ...
    'lambda', 0.002 );

material = MaterialClass.preset(3);

observation = struct('x', 0:0.1:20, ... % size of bins in space
    'y', [-pi pi], ...
    'z', [-pi/2 pi/2], ...
    'directions', [0 pi], ...
    'time', 0:0.01:10 );

d = geometry.dimension;
vp = material.vp; vs = material.vs;
K = vp/vs;
Sigmapp = 0.2*vp; Sigmaps = 0.2*vp; Sigmap = Sigmapp + Sigmaps;
Sigmasp = Sigmaps/((d-1)*K^d); Sigmass = Sigmap -Sigmasp;

material.sigma = {@(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmapp, @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmaps; ...
    @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmasp, @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmass};

Wsp = 0; % S to P wave energy ratio at the source

inds = [20 50 80]; % index of the desired observation points

% running our code, Monte Carlo-based
obs = radiativeTransfer( geometry, source, material, observation );

Ep = squeeze(obs.energyDensity(:,:,:,1))/obs.dz;
Es = squeeze(obs.energyDensity(:,:,:,2))/obs.dz;
Etotus = Ep + Es;

% running Yoshimoto's Monte Carlo-based approach
EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );
EtotY = sum(EY,3);

% computing semi-analytical solution (Nakahara 2011, Sato 1994)
Eanalytical = Comparison.analyticalEnergyIsotropicElastic(geometry, material, observation, obs.x(inds), Wsp);

% comparison of P energy densities
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Ep(inds,:), '-k' );
h2 = semilogy( obs.t, EY(inds,:,1), '-b' );
set(gca, 'YScale', 'log');
legend( [h1(1), h2(1)], {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('P-wave energy densities at different source-station distances')
title(titlecase);

% comparison of S energy densities
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Es(inds,:), '-k' );
h2 = semilogy( obs.t, EY(inds,:,2), '-b' );
set(gca, 'YScale', 'log');
legend( [h1(1), h2(1)], {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('S-wave energy densities at different source-station distances')
title(titlecase);

% comparison of total energy densities
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Etotus(inds,:), '-k' );
h2 = semilogy( obs.t, EtotY(inds,:), '-r' );
h3 = semilogy( obs.t, Eanalytical, '-b' );
set(gca, 'YScale', 'log');
legend( [h1(1), h2(1), h3(1)], {'Monte Carlo (our code)', ...
    'Monte Carlo (Yoshimoto 2000)','Analytical (Sato, 1994)'},'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Energy densities at different source-station distances')
currentYLim = ylim(gca);
ylim(gca, [1e-5 currentYLim(2)]);
title(titlecase);