% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

% Physics
physics = struct( 'acoustics', true, ...
                  'dimension', 2 );
                   
% Point source
source = struct( 'numberParticles', 1e7, ...
                 'position', [0 0 -2], ... 
                 'lambda', 0.1 );
             
% material properties
material = struct( 'v', 1, ...
                   'sigma', @(th) 1/4/pi*ones(size(th)));
               
% observations
observation = struct('dx', 0.05, ...           % size of bins in space
                     'time', 0:0.1:100, ...     % observation times
                     'Ndir', 100 );            % number of bins for directions

% format of 'plane' arguments is [ z type ] (possibly two lines)
%           z indicates the altitude of the plane along the axis
%           type is -1 for Dirichlet and 1 for Neumann 
% the movieBox is always in [Lx Lz], y=0. When a plane is defined, the box
% sticks to that plane, otherwise, it is [-L/2 L/2].
geometry = struct( 'planeX', [-5 1;1.5 1], ...
                   'planeZ', [0 1;-3 1], ...
                   'movieBox', [10 7]);

% radiative transfer solution - 2D - acoustic
obs = radiativeTransfer( physics, source, material, observation, geometry );
%M = plotGrid('full',obs,1);
M = plotGrid('half',obs,0.5);
%scatterDirections(obs,30);
