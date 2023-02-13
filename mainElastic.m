% This code works in 2D/3D and assumes the source is centered on 0. The
% initial distance from 0 is modeled gaussian with standard deviation 
% source.lambda and uniformly-distributed angle(s). The initial direction is 
% uniform.

% Physics
physics = struct( 'acoustics', true, ...
                  'dimension', 2 );

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'lambda', 0.1, ...
                 'polarization', 'P');  % 'P' or 'S'
             
% material properties
material = struct( 'acoustics', false, ...
                   'dimension', 2, ...
                   'vp', 1.5, ...
                   'vs', 1, ...
                   'sigmaPP', @(th) 1/8/pi*ones(size(th)), ...
                   'sigmaPS', @(th) 2*1.5^2/20/pi*ones(size(th)), ...
                   'sigmaSP', @(th) 1/20/pi*ones(size(th)), ...
                   'sigmaSS', @(th) 1/8/pi*ones(size(th)) );
               
% observations
observation = struct('sensors', 0:0.05:10, ... % bins for histograms
                     'time', 0:0.1:10, ...     % observation times
                     'Ndir', 100 );            % direction discretization

% radiative transfer solution - 2D - elastic
%obs = radiativeTransferElastic( source, material, observation );
obs = radiativeTransfer( physics, source, material, observation );
M = plotGrid('full',obs,1);
scatterDirections(obs,30);
