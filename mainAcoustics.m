% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

% geometry
geometry = struct( 'dimension', 3 , ...
                   'frame', 'cartesian' ); % 'spherical' (default) or 'cartesian'
% add boundaries, specifying fields 'dir' for normal direction to the plane
% (1='x', 2='y', 3='z'), and 'val' for the position of the plane
geometry.bnd(1) = struct('dir',3,'val',0);

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'type', 'plane', ...         % 'point' (default) or 'plane'
                 'position', [0 0 -2], ...    % always in cartesian frame
                 'direction', 3,          ... % only used with point sources ('uniform' (default) or 'outgoing')
                 'radial', 0, ...             % only used with point sources. Initial condition function of radius (always spherical coordinates)
                 'extent', [10 10], ...       % only used with plane sources. Extent of the source around 'position' (always cartesian coordinates)
                 'lambda', 0.1 );    

% observations
% choose 2 variables only to perform histograms (among x, y, z, directions)
% if geometry.frame = 'spherical', (x,y,z) correspond respectively to
%  x=r, y=azimuth (in [-pi pi]), z=elevation (in [-pi/2 pi/2])
observation = struct('x', -2:.1:2, ...                  % bins in space
                     'y', [-Inf Inf], ...                 
                     'z', -4:0.1:0, ...                 % unused in 2D
                     'directions', [0 pi], ...          % bins for directions [0 pi]         
                     'time', 0:0.05:3 );                % observation times

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

plotting = struct( 'checkEnergy', false, ...
                   'movieTotalEnergy', true, ...
                   'timehistory', false, ...
                   'sensors', sensors );

plotEnergies( plotting, obs, material, source.lambda )