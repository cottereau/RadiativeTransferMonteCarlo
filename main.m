% This code works in 2D and assumes the source is centered on [0 0]. The
% initial distance from (0,0) is modeled gaussian with standard deviation 
% source.lambda and uniformly-distributed angle. The initial direction is 
% uniform.

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'lambda', 0.1);
             
% material properties
material = struct( 'acoustics', true, ...
                   'dimension', 2, ...
                   'v', 1, ...
                   'sigma', @(th) 1/4/pi*ones(size(th)));
               
% observations
observation = struct('sensors', 0.001:0.05:10, ... % to observe at given sensors
                     'check', false, ... % to perform some checks for debug
                     'time', 0:0.1:10, ... % times at which to evaluate passage in points
                     'Ndir', 100, ... % number of directions to be considered
                     'movie', true ); % create a movie from the sensors (using rotational invariance)

% geometry of the half-space problem (defined by z<0)
geometry = struct('sourcePosition', -2, ... % vertical position of the source
                  'boundaryCondition', 'Neumann' ); % 'Neumann' or 'Dirichlet'

% radiative transfer solution - 2D - acoustic
%obs = radiativeTransferHS( source, material, observation, geometry );
obs = radiativeTransfer( source, material, observation );
%plotPoints(obs)
   M = plotGrid('full',obs,0.05);
%M = plotGrid('half',obs,0.01);
scatterDirections(obs,100);
