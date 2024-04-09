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
                   'dimension', 3 );

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'position', [2.5 1 -2], ... 
                 'lambda', 0.001 );

% observations
observation = struct('dr', 0.05, ...        % size of bins in space
                     'time', 0:0.05:20, ...  % observation times
                     'Ndir', 10 );         % number of bins for directions           

% material properties
% material.coefficients_of_variation defines the coefficients of variation
% of kappa (bulk modulus) and rho (density), respectively.
% material.correlation_coefficients defines the correlation coefficient
% between kappa (bulk modulus) and rho (density).
material = struct( 'acoustics', true, ...
                   'v', 1, ...
                   'Frequency', 10, ...
                   'correlationLength', 10, ...
                   'spectralType', 'exp', ...
                   'coefficients_of_variation', [0.1 0.2], ...
                   'correlation_coefficients', -0.5 );

% radiative transfer solution - acoustic with boundaries
switch geometry.type
    case 'fullspace'
        obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
    otherwise
        obs = radiativeTransferAcoustics( source, material, observation, geometry );
end

% plotting output
sensors = [3   1 -0.5; 
           3.9 1 -2.5;
           0.2 1 -1.5;
           1.5 1 -2.5;
           3.5 1 -2.8];

plotting = struct( 'equipartition', true, ...
                   'movieTotalEnergy', true, ...
                   'movieDirectionalEnergy', true, ...
                   'sensors', sensors);

plotEnergies( plotting, obs, material, source.lambda )