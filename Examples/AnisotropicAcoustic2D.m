close all
clear all
clc

titlecase = '2D acoustic case with anisotropic scattering';
disp(['Testing ' titlecase ' ...']);

% input data
geometry = struct( 'dimension', 2 );

source = struct( 'numberParticles', 1e6, ...
    'position', [0 0], ...
    'lambda', 2e-4);

material = MaterialClass.preset(1);
% forward scattering regime
material.sigma = {@(th) 1/4/pi*(1+4*cos(th).^4)};

observation = struct( 'x', 0:0.03:9, ...
    'y', [-pi pi], ...
    'directions', [0 pi], ...
    'time', 0:0.02:20 );

inds = [60 150 240]; % index of the desired observation points

% running our code, Monte Carlo-based
obs = radiativeTransferUnbounded( geometry, source, material, observation );
Eus = squeeze(obs.energyDensity);

% running Yoshimoto's Monte Carlo-based approach
EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );

% visual comparison
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Eus(inds,:), '-k' );
h2 = semilogy( observation.time, EY(inds,:), '-r' );
set(gca, 'YScale', 'log'); ylim([1e-5 1]);
legend( [h1(1), h2(1)], ...
    {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);