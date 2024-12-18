% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

% geometry
geometry = struct( 'dimension', 3 , ...
                   'frame', 'cartesian' ); % 'spherical' (default) or 'cartesian'

% Point source
r = [0 5 7 9 11 13 15 17 19 10000]/100;
s = [0 .0743 .9697 .2064 .0571 .0544 .0013 .0066 0 0]*1e-13;
s = griddedInterpolant(r,s);
source = struct( 'numberParticles', 1e6, ...
                 'position', [2.5 1 -2], ...  % always in cartesian frame
                 'lambda', 0.1, ...
                 'direction', 'outgoing', ... % initial direction 'uniform' (default) or 'outgoing' or 'plane'
                 'radial', s );               % initial condition function of radius (always spherical)

% observations
% choose 2 variables only to perform histograms (among x, y, z, directions)
% if geometry.frame = 'spherical', (x,y,z) correspond respectively to
%  x=r, y=azimuth (in [-pi pi]), z=elevation (in [-pi/2 pi/2])
observation = struct('x', 0:.1:10, ...                    % bins in space
                     'y', [-Inf Inf], ...                 
                     'z', [-Inf Inf], ...                 % unused in 2D
                     'directions', linspace(0,pi,10), ... % number of bins for directions [0 pi]         
                     'time', 0:0.3:1 );                   % observation times

% material properties - more examples can be found in Example folder
% here is a basic example
freq = 10; % in Hz
material = MaterialClass( geometry, ...
                          freq, ...
                          true, ...          % true for acoustics
                          1, ...             % average wave velocity
                          [0.1 0.2], ...     % coefficients of variation of kappa and rho.
                          -0.5, ...          % correlation coefficient of kappa/rho
                          'exp', ...         % autocorrelation function
                           0.1 );            % correlation length
material.timeSteps = 0;                      % time Steps : 0=small 1=large
                          
% radiative transfer solution - acoustic with boundaries
obs = radiativeTransferUnbounded( geometry, source, material, observation );

% plotting output
sensors = [0 0 0; 
           0 0 5];

plotting = struct( 'checkEnergy', true, ...
                   'movieEnergy', false, ...
                   'timehistory', false, ...
                   'sensors', sensors);

plotEnergies( plotting, obs, material, source.lambda )