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
observation = struct('sensors', 0:0.05:10, ... % bins for histograms
                     'time', 0:0.1:10, ...     % observation times
                     'Ndir', 100 );            % direction discretization

% geometry of the half-space problem (defined by z<0)
geometry = struct('sourcePosition', -2, ... % vertical position of the source
                  'boundaryCondition', 'Neumann' ); % 'Neumann' or 'Dirichlet'

% radiative transfer solution - 2D - acoustic
%obs = radiativeTransferHS( source, material, observation, geometry );
obs = radiativeTransfer( source, material, observation );
M = plotGrid('full',obs,0.05);
%M = plotGrid('half',obs,0.01);
scatterDirections(obs,30);
