close all
clear all
clc

titlecase = '2D acoustic case with isotropic scattering';
disp(['Testing ' titlecase ' ...']);

% input data
geometry = struct( 'dimension', 2 );

source = struct( 'numberParticles', 1e6, ...
    'position', [0 0], ...
    'lambda', 2e-4 );
material = MaterialClass.preset(1);
observation = struct('x', 0:0.03:9, ...
    'y', [-pi pi], ...
    'directions', [0 pi], ...
    'time', 0:0.05:20 );

inds = [60 150 240]; % index of the desired observation points

% running our code, Monte Carlo-based
obs = radiativeTransferUnbounded( geometry, source, material, observation );
Eus = squeeze(obs.energyDensity);

% running Yoshimoto's Monte Carlo-based approach
EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );

% computing Paasschens solution
[EP,Ediff] = Comparison.analyticalPaasschens( material, observation, geometry );

% visual comparison
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Eus(inds,:), '-k' );
h2 = semilogy( obs.t, EY(inds,:), '-r' );
h3 = semilogy( obs.t, EP(inds,:), '-b' );
h4 = semilogy( obs.t, Ediff(inds,:), ':r' );
set(gca, 'YScale', 'log'); ylim([1e-5 1]);
legend( [h1(1), h2(1), h3(1), h4(1)], ...
    {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
    'Analytical (Paasschens, 1997)', 'Diffusion approximation'}, ...
    'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);