% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

geometry = struct( 'dimension', 2 );

% Point source
source = struct( 'numberParticles', 1e7, ...
                 'polarization', 'P', ...
                 'lambda', 0.1 );

% observations
observation = struct('dr', 0.05, ...        % size of bins in space
                     'time', 0:5:100, ...  % observation times
                     'Ndir', 30 );         % number of bins for directions           
 
% material properties
% material.coefficients_of_variation defines the coefficients of variaiton
% of lambda, mu (Lamé coefficients) and rho (density), respectively.
% material.correlation_coefficients defines the correlation coefficient
% between (lambda,mu), (lambda,rho), and (mu,rho), respectively.
material = struct( 'acoustics', false, ...
                   'vp', sqrt(2), ...
                   'vs', 1, ...
                   'Frequency', 2*pi/source.lambda, ...
                   'correlationLength', 0.2, ...
                   'spectralType', 'exp', ...
                   'coefficients_of_variation', [0.1 0.2 0.3], ...
                   'correlation_coefficients', [-0.1 0.15 -0.5]);
material.sigma = PSDF2sigma( geometry.dimension, material );

% radiative transfer solution - acoustic with boundaries
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% plotting output
sensors = [1 0 0];
plotEnergies( obs, material, source.lambda, sensors );
