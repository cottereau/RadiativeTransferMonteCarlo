close all
clear all
clc

titlecase = '3D acoustic case with anisotropic scattering';
disp(['Testing ' titlecase ' ...']);

% input data
geometry = struct( 'dimension', 3 );

source = struct( 'numberParticles', 1e6, ...
    'position', [0 0 0], ...
    'lambda', 2e-4 );

freq = 10; % in Hz
material = MaterialClass( geometry, ...
    freq, ...
    true, ...          % true for acoustics
    1, ...             % average wave velocity
    [0.1 0.2], ...     % coefficients of variation of kappa and rho.
    -0.5, ...          % correlation coefficient of kappa/rho
    'exp', ...         % autocorrelation function
    0.1);              % correlation length

observation = struct( 'x', 0:0.1:10, ...
    'y', [-pi pi], ...
    'z', [-pi/2 pi/2], ...
    'directions', [0 pi], ...
    'time', 0:0.02:20, ...  % observation times
    'Ndir', 10 );           % number of bins for directions

inds = [20 50 80]; % index of the desired observation points

% running our code, Monte Carlo-based
obs = radiativeTransfer( geometry, source, material, observation );
Eus = squeeze(obs.energyDensity)/obs.dz;

% running Yoshimoto's Monte Carlo-based approach
EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );

% comparison
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Eus(inds,:), '-k' );
h2 = semilogy( obs.t, EY(inds,:), '-r' );
set(gca, 'YScale', 'log'); ylim([1e-5 1]);
legend( [h1(1), h2(1)], ...
    {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);