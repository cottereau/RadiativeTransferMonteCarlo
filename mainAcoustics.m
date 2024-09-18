% Approximate solution of the Radiative Transfer Equation using the 
% Monte Carlo method

% geometry
geometry = struct( 'dimension', 3 );

% source
r = [0 5 7 9 11 13 15 17 19 10000];
s = [0 .0743 .9697 .2064 .0571 .0544 .0013 .0066 0 0]*1e-13;
s = griddedInterpolant(r,s);
source = struct( 'numberParticles', 1e7, ...
                 'position', [2.5 1 -2], ... 
                 'radial', s );      % initial condition function of radius

% material properties - more examples can be found in Example folder
% here is a basic example
freq = 1; % in rad/s
material = MaterialClass( geometry, ...
                          freq, ...
                          true, ...          % true for acoustics
                          1, ...             % average wave velocity
                          [0.1 0.2], ...     % coefficients of variation of kappa and rho.
                          -0.5, ...          % correlation coefficient of kappa/rho
                          'exp', ...         % autocorrelation function
                          0.1);              % correlation coefficient

% observations
% choose two variables only among (r, theta, phi, psi) to perform histograms 
observation = struct('r', 0:.1:1, ...                 % bins in space
                     'azimuth', linspace(-pi,pi,10), ... % bins for azimuth angle [-pi pi]
                     'elevation', [-pi/20 pi/20], ... % bins for elevation angle [-pi/2 pi/2] in 3d only
                     'directions', [0 pi], ...        % number of bins for directions [0 pi]         
                     'time', 0:0.3:10 );              % observation times
                     
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