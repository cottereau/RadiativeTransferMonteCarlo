close all
clearvars
clc

titlecase = 'Movie for the 3D elastic case with anisotropic scattering';
disp(['Testing ' titlecase ' ...']);

% Input geometry
geometry = struct( ...
    'dimension', 3, ...
    'frame', 'cylindrical');

% Point source initially emitting P waves
source = struct( ...
    'numberParticles', 5e6, ...
    'position', [0 0 0], ...
    'polarization', 'P', ...
    'lambda', 0.002);

% The following setup favors a stochastic scattering regime
freq = 10; % in Hz
material = MaterialClass(geometry, ...
    freq, ...
    false, ...            % false for elasticity
    [6 6/sqrt(3)], ...    % P- and S-wave velocities
    [0.1 0.1 0.], ...     % coefficients of variation of lambda, mu and rho
    [0.1 0. 0.], ...      % correlations: lambda/mu, lambda/rho and mu/rho
    'exp', ...             % autocorrelation function
    0.1);                 % correlation length

% No intrinsic attenuation
material.Q = [Inf Inf];

% Observe energy in the cylindrical r-z plane
observation = struct( ...
    'x', 0:0.2:20, ...         % cylindrical radius r
    'y', [-pi pi], ...         % integrate over azimuth
    'z', -20:0.2:20, ...       % resolve z
    'directions', [0 pi], ...  % integrate over propagation direction
    'time', 0:0.02:3);

% Run the Monte Carlo radiative-transfer simulation
obs = radiativeTransfer(geometry, source, material, observation);

% Generate connected P-, S- and total-energy panels
plotting = struct( ...
    'movieTotalEnergy', true, ...
    'energyScale', 'linear', ... % or log10
    'filename', 'AnisotropicElastic3D_rz.gif', ...
    'delayTime', 0.05, ...
    'visible', 'on', ...
    'axisMode', 'equal');

plotEnergies(plotting, obs, material, source.lambda);