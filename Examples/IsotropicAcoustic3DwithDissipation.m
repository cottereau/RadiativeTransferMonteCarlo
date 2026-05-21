close all
clearvars
clc

titlecase = '3D acoustic isotropic scattering with Q-based absorption';
disp(['Testing ' titlecase ' ...']);

% Input data
geometry = struct('dimension', 3);

source = struct('numberParticles', 1e6, ...
                'position', [0 0 0], ...
                'lambda', 2e-4);

material = MaterialClass.preset(1);
material.Frequency = 10; % Hz

% Quality factor
material.Q = 100;

observation = struct('x', 0:0.1:10, ...
                     'y', [-pi pi], ...
                     'z', [-pi/2 pi/2], ...
                     'directions', [0 pi], ...
                     'time', 0:0.05:20);

inds = [20 50 80];

% running our code, Monte Carlo-based
obs = radiativeTransfer(geometry, source, material, observation);

Eus = squeeze(obs.energyDensity);

% computing Paasschens solution
[EP, Ediff] = Comparison.analyticalPaasschens(material, observation, geometry);

% visual comparison
figure; hold on; grid on; box on;

h1 = semilogy(obs.t, Eus(inds,:), '-k', 'LineWidth', 1.1);
h2 = semilogy(obs.t, EP(inds,:), '-b', 'LineWidth', 1.1);
h3 = semilogy(obs.t, Ediff(inds,:), ':r', 'LineWidth', 1.2);

set(gca, 'YScale', 'log');
ylim([1e-8 1]);

legend([h1(1), h2(1), h3(1)], ...
    {'Monte Carlo (our code)', ...
     'Analytical Paasschens with dissipation', ...
     'Diffusion approximation with dissipation'}, ...
     'FontSize', 12, 'Location', 'best');

xlabel('Lapse time [s]');
ylabel('Integrated energy density');
title(titlecase);