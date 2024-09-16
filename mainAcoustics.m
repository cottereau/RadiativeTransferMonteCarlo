% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

% geometry
geometry = struct( 'dimension', 3 );

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'position', [2.5 1 -2], ... 
                 'lambda', 0.1 ); %wavelength 

% observations
% choose two variables only to perform histograms (r, theta, phi, psi)
observation = struct('r', 0:.1:1, ...                 % bins in space
                     'azimuth', linspace(-pi,pi,10), ... % bins for azimuth angle [-pi pi]
                     'elevation', [-pi/20 pi/20], ... % bins for elevation angle [-pi/2 pi/2] in 3d only
                     'directions', [0 pi], ...        % number of bins for directions [0 pi]         
                     'time', 0:0.3:10 );              % observation times

% material properties - more examples can be found in Example folder
% here is a basic example
material = MaterialClass( geometry, ...
                          source, ...
                          true, ...          % true for acoustics
                          1, ...             % average wave velocity
                          [0.1 0.2], ...     % coefficients of variation of kappa and rho.
                          -0.5, ...          % correlation coefficient of kappa/rho
                          'exp', ...         % autocorrelation function
                          0.1);              % correlation coefficient
                          
% radiative transfer solution - acoustic with boundaries
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% plotting output
sensors = [0 0 0; 
           0 0 5];

plotting = struct( 'equipartition', false, ...
                   'movieTotalEnergy', false, ...
                   'movieDirectionalEnergy', false, ...
                   'timehistory', true, ...
                   'sensors', sensors);

plotEnergies( plotting, obs, material, source.lambda )