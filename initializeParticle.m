function P = initializeParticle( N, d, acoustics, source )

% initial position of each particle: radius follows a Gaussian law with
% standard deviation lambda, and angle follows a uniform law.
r = abs(randn(N,1)*source.lambda);

% initial angles and directions
if d==2
    phi0 = (pi/2)*ones(N,1);
    phi = (pi/2)*ones(N,1);
elseif d==3
    phi0 = 2*pi*rand(N,1);
    phi = 2*pi*rand(N,1);
end
theta0 = 2*pi*rand(N,1);
theta = 2*pi*rand(N,1);
dir = [cos(theta).*sin(phi) sin(theta).*sin(phi) cos(phi)];
perp = [-sin(theta) cos(theta) zeros(N,1)];

% initial positions in cartesian coordinates
x = r.*[cos(theta0).*sin(phi0) sin(theta0).*sin(phi0) cos(phi0)];

% current time for the particle
t = zeros(N,1);

% initial polarisation of each particle
% in acoustics, this variable is unused (and set always to true)
% in elastics, true corresponds to P waves (default), and false to S waves
p = true(N,1);
if ~acoustics && isfield(source,'polarization') && source.polarization=='S'
    p = false(N,1);
end

% coherent flag (false when particle has been scattered at least once)
coherent = true(N,1);
 
% initialize structure
P = struct( 'd', d, ...                 % dimension of the problem
            'N', N, ...                 % number of particles
            'x', x, ...                 % cartesian coordinates
            'dir', dir, ...             % direction of propagation
            'perp', perp, ...           % orthogonal to direction of propagation
            'p', p, ...                 % polarization (used only in elasticity)
            't', t, ...                 % current time for the particle
            'coherent', coherent );     % false when particle has been scattered

