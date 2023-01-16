% this code works in 2D and assumes invariance by rotation around the 
% source, so the source is always assumed centered on [0 0]. Its shape is a
% Gaussian (cylindrical) with typical size lambda, and can have any 
% polarization.
%
% Point source
source = struct( 'numberParticles', 100000, ...
                 'lambda', 0.1, ...
                 'polarization', 'P');  % 'P' or 'S'
             
% material properties
material = struct( 'vp', 1.5, ...
                   'vs', 1, ...
                   'sigmaPP', @(th) 1/8/pi*ones(size(th)), ...
                   'sigmaPS', @(th) 2*1.5^2/20/pi*ones(size(th)), ...
                   'sigmaSP', @(th) 1/20/pi*ones(size(th)), ...
                   'sigmaSS', @(th) 1/8/pi*ones(size(th)) );
               
% observations
observation = struct('sensors', 0.001:0.05:10, ... % to observe at given sensors
                     'check', false, ... % to perform some checks for debug
                     'time', 0:0.1:10, ... % times at which to evaluate passage in points
                     'movie', true ); % create a movie from the sensors (using rotational invariance)

% radiative transfer solution - 2D - elastic
obs = radiativeTransferElastic( source, material, observation );
plotPoints(obs)
M = plotGrid('full',obs,0.01);

