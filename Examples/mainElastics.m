% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.
clc

geometry = struct( 'dimension', 3 );

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'polarization', 'P', ...
                 'lambda', 0.1 );

% observations
% choose 2 variables only to perform histograms (among x, y, z, directions)
% if geometry.frame = 'spherical', (x,y,z) correspond respectively to
%  x=r, y=azimuth (in [-pi pi]), z=elevation (in [-pi/2 pi/2])
observation = struct('x', 0:.04:5, ...                    % bins in space
                     'y', [-Inf Inf], ...                 
                     'z', [-Inf Inf], ...                 % unused in 2D
                     'directions', linspace(0,pi,10), ... % number of bins for directions [0 pi]         
                     'time', 0:0.1:1 );                   % observation times

% material properties - more examples can be found in Example folder
% here is a basic example
freq = 10; % in Hz
material = MaterialClass( geometry, ...
                          freq, ...
                          false, ...            % true for acoustics
                          [6 6/sqrt(3)], ...    % defines the velocity of pressure waves and the shear waves
                          [0.8 0.8 0.], ...     % defines the coefficients of variation of lambda, mu (Lamé coefficients) and rho (density), respectively.
                          [0.1 0. 0.], ...      % defines the correlation coefficient between (lambda,mu), (lambda,rho), and (mu,rho), respectively
                          'exp', ...            % defines the autocorrelation function
                          0.1);                 % defines the correlation length
material.timeSteps = 0;                         % time Steps : 0=small 1=large

% radiative transfer solution - acoustic with boundaries
obs = radiativeTransfer( geometry, source, material, observation );

% plotting output
plotting = struct( 'equipartition', true, ...
                   'movieTotalEnergy', true, ...
                   'movieDirectionalEnergy', false, ...
                   'timehistory', true, ...
                   'sensors', [1 0 0; 9 0 0; 9 0 1; 9 1 0]);

plotEnergies( plotting, obs, material, source.lambda );
