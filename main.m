% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.

% - 'type' is either 'fullspace', 'halfspace', 'slab', or 'box'.
% 'halfspace' is defined for z<0. 'slab' is bounded between two planes of
% constant z. In 2D, 'box' is bounded in x and z.
% - 'size' is a vector of length the dimension of the problem. It specifies
% the [width depth height] of the simulation and/or plotting box. More
% specifically:
%  (1) height is the size (in z) of the plotting box, and also the size of
%      the simulation box, in the case of types 'slab' and 'box'. The
%      range of coordinates of the box in z are [-heigth 0]
%  (2) width is the size (in x) of the plotting box, and also the length
%      of the simulation box, in the case of type 'box'. The range of
%      coordinates of the box in x are [0 width]
%  (3) depth is the size (in y) of the simulation box for type 'box' in 3
%      dimensions (unused for other cases).  The range of coordinates of
%      the box in y are [0 depth]
% The plotting box is always drawn in the plane of the source (in y). For
% now, only homogeneous Neumann boundary conditions are enforced
geometry = struct( 'type', 'box', ...
                   'size', [4 3 3], ...
                   'dimension', 2 );

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'position', [2.5 1 -2], ... 
                 'lambda', 0.1 );

% observations
observation = struct('dr', 0.05, ...        % size of bins in space
                     'time', 0:0.05:2, ...  % observation times
                     'Ndir', 5 );         % number of bins for directions           
 
% material properties
material = struct( 'acoustics', true, ...
                   'v', 1, ...
                   'Frequency', 2*pi/source.lambda, ...
                   'correlationLength', 10, ...
                   'spectralType', 'exp', ...
                   'correlationMatrix', 0.1*ones(2) );
material.sigma = PSDF2sigma( geometry.dimension, material );

% radiative transfer solution - acoustic with boundaries
% obs = radiativeTransferAcoustics( source, material, observation, geometry );
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% plotting output
sensors = [3   1 -0.5; 
           3.9 1 -2.5;
           0.2 1 -1.5;
           1.5 1 -2.5;
           3.5 1 -2.8];
plotEnergies( obs, material, source.lambda, 4, sensors );
