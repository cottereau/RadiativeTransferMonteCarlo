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
observation = struct('dr', 0.04, ...        % size of bins in space
                     'time', 0:0.1:15, ...  % observation times
                     'Ndir', 10 );         % number of bins for directions           

% material properties - more examples can be found in Example folder
% here is a basic example
material = MaterialClass( geometry, ...
                          source, ...
                          false, ...            % true for acoustics
                          [6 6/sqrt(3)], ...    % defines the velocity of pressure waves and the shear waves
                          [0.8 0.8 0.], ...     % defines the coefficients of variation of lambda, mu (Lamé coefficients) and rho (density), respectively.
                          [0.1 0. 0.], ...      % defines the correlation coefficient between (lambda,mu), (lambda,rho), and (mu,rho), respectively
                          'exp', ...            % defines the autocorrelation function
                          0.1);                 % defines the correlation coefficient

% radiative transfer solution - acoustic with boundaries
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% plotting output
plotting = struct( 'equipartition', true, ...
                   'movieTotalEnergy', false, ...
                   'movieDirectionalEnergy', false, ...
                   'timehistory', true, ...
                   'sensors', [1 0 0; 9 0 0; 9 0 1; 9 1 0]);

plotEnergies( plotting, obs, material, source.lambda );
