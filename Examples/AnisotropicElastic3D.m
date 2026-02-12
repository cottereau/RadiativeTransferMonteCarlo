close all
clear all
clc

titlecase = '3D elastic case with anisotropic scattering';
disp(['Testing ' titlecase ' ...']);

geometry = struct( 'dimension', 3 );

source = struct( 'numberParticles', 1e6, ...
    'position', [0 0 0], ...
    'polarization', 'P', ...
    'lambda', 0.002 );

% The following setup favors a stochastic scattering regime
freq = 10; % in Hz
material = MaterialClass( geometry, ...
    freq, ...
    false, ...            % true for acoustics
    [6 6/sqrt(3)], ...    % defines the velocity of pressure waves and the shear waves
    [0.8 0.8 0.], ...     % defines the coefficients of variation of lambda, mu (Lamé coefficients) and rho (density), respectively.
    [0.1 0. 0.], ...      % defines the correlation coefficient between (lambda,mu), (lambda,rho), and (mu,rho), respectively
    'exp', ...            % defines the autocorrelation function
    0.1);                 % defines the correlation length
material = prepareSigma( material, geometry.dimension );

observation = struct('x', 0:0.1:20, ... % size of bins in space
    'y', [-pi pi], ...
    'z', [-pi/2 pi/2], ...
    'directions', [0 pi], ...
    'time', 0:0.01:10 );

inds = [20 50 80]; % index of the desired observation points

% running our code, Monte Carlo-based
obs = radiativeTransfer( geometry, source, material, observation );

Ep = squeeze(obs.energyDensity(:,:,:,1))/obs.dz;
Es = squeeze(obs.energyDensity(:,:,:,2))/obs.dz;
Etotus = Ep + Es;

% running Yoshimoto's Monte Carlo-based approach
EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation, false );
EtotY = sum(EY,3);

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
set(gca, 'YScale', 'log');
legend( [h1(1), h2(1)], {'Monte Carlo (our code)','Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Energy densities at different source-station distances')
currentYLim = ylim(gca);
ylim(gca, [1e-5 currentYLim(2)]);
title(titlecase);