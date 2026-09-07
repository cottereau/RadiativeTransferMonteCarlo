close all
clearvars
clc

titlecase = 'Movie for the 3D acoustic case with anisotropic scattering';
disp(['Testing ' titlecase ' ...']);

% input data
geometry = struct( ...
    'dimension', 3, ...
    'frame', 'cylindrical');

source = struct( 'numberParticles', 5e6, ...
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

% No intrinsic attenuation
material.Q = Inf;

observation = struct( ...
    'x', 0:0.05:6, ...          % cylindrical radius r
    'y', [-pi pi], ...         % integrate over azimuth
    'z', -6:0.05:6, ...         % resolve z
    'directions', [0 pi], ...  % integrate over propagation direction
    'time', 0:0.05:5);

% running our code, Monte Carlo-based
obs = radiativeTransfer( geometry, source, material, observation );

plotting = struct( ...
    'movieTotalEnergy', true, ...
    'energyScale', 'linear', ... % or log10
    'filename', 'AnisotropicAcoustic3D_rz.gif', ...
    'delayTime', 0.05, ...
    'visible', 'on', ...
    'axisMode', 'equal');

plotEnergies(plotting, obs, material, source.lambda);