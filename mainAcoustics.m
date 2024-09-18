% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

% geometry
geometry = struct( 'dimension', 3 );

% Point source
r = [0 5 7 9 11 13 15 17 19 10000]/100;
s = [0 .0743 .9697 .2064 .0571 .0544 .0013 .0066 0 0]*1e-13;
s = griddedInterpolant(r,s);
source = struct( 'numberParticles', 1e6, ...
                 'position', [2.5 1 -2], ... 
                 'lambda', 0.1, ...
                 'radial', s );      % initial condition function of radius

% observations
% choose two variables only to perform histograms (r, theta, phi, psi)
observation = struct('r', 0:.1:10, ...                 % bins in space
                     'azimuth', linspace(-pi,pi,10), ... % bins for azimuth angle [-pi pi]
                     'elevation', [-pi/20 pi/20], ... % bins for elevation angle [-pi/2 pi/2] in 3d only
                     'directions', [0 pi], ...        % number of bins for directions [0 pi]         
                     'time', 0:0.3:1 );              % observation times

% material properties - more examples can be found in Example folder
% here is a basic example
freq = 2*pi*10; % in rad/s
material = MaterialClass( geometry, ...
                          freq, ...
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

plotting = struct( 'checkEnergy', true, ...
                   'movieEnergy', false, ...
                   'timehistory', false, ...
                   'sensors', sensors);

plotEnergies( plotting, obs, material, source.lambda )