% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

geometry = struct( 'dimension', 2 );

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'polarization', 'P', ...
                 'lambda', 0.1 );

% observations
observation = struct('dr', 0.05, ...        % size of bins in space
                     'time', 0:0.05:2, ...  % observation times
                     'Ndir', 5 );         % number of bins for directions           
 
% material properties
material = struct( 'acoustics', false, ...
                   'vp', sqrt(2), ...
                   'vs', 1, ...
                   'Frequency', 2*pi/source.lambda, ...
                   'correlationLength', 10, ...
                   'spectralType', 'exp', ...
                   'correlationMatrix', 0.1*eye(3) );
material.sigma = PSDF2sigma( geometry.dimension, material );

% radiative transfer solution - acoustic with boundaries
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% plotting output
sensors = [3   1 -0.5; 
           3.9 1 -2.5;
           0.2 1 -1.5;
           1.5 1 -2.5;
           3.5 1 -2.8];
plotEnergies( obs, material, source.lambda, 4, sensors );
