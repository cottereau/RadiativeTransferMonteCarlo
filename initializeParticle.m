function P = initializeParticle( N, d, acoustics, source, material )

% initial position of each particle: radius follows a Gaussian law with
% standard deviation lambda, and angle follows a uniform law.
r = abs(randn(N,1)*source.lambda);

% initial angles and directions
if d==2
    phi = (pi/2)*ones(N,1);
elseif d==3
    phi = 2*pi*rand(N,1);
end
theta = 2*pi*rand(N,1);
dir = [cos(theta).*sin(phi) sin(theta).*sin(phi) cos(phi)];
perp = [-sin(theta) cos(theta) zeros(N,1)];

% initial positions in cartesian coordinates
x = [r.*cos(theta).*sin(phi) r.*sin(theta).*sin(phi) r.*cos(phi)];

% initial polarisation of each particle
% in acoustics, this variable is unused (and set always to true)
% in elastics, true corresponds to P waves (default), and false to S waves
p = true(N,1);
if isfield(source,'polarization') && source.polarization=='S'
    p = false(N,1);
end

% current time for the particle
t = zeros(N,1);

% material velocity of the background for each particle
if acoustics && isfield(material,'v')
    v = material.v*ones(N,1);
elseif ~physics.acoustics && isfield(material,'vp') && isfield(material,'vs')
    v = material.vs*ones(N,1);
    v(p) = material.vp;
else
    error('unknown velocity of the background')
end

% meanFreePath for each particle
if isfield(material,'meanFreePath')
    mfp = material.meanFreePath*ones(N,1);
elseif isfield(material,'meanFreePathP') && isfield(material,'meanFreePathS')
    mfp = material.meanFreePathS*ones(N,1);
    mfp(p) = material.meanFreePathP;
else
    error('unknown meanFreePath')
end

% initialize structure
P = struct( 'd', d, ...                 % dimension of the problem
            'acoustics', acoustics, ... % true=acoustics, false=elastics
            'N', N, ...                 % number of particles
            'x', x, ...                 % cartesian coordinates
            'dir', dir, ...             % direction of propagation
            'perp', perp, ...           % orthogonal to direction of propagation
            'p', p, ...                 % polarization (used only in elasticity)
            'meanFreePath', mfp, ...    % mean free path
            'v', v, ...                 % propagation velocity
            't', t );                   % current time for the particle

