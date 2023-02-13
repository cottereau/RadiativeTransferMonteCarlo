function P = initializeParticle( N, physics, source, material )

% initial position of each particle: radius follows a Gaussian law with
% standard deviation lambda, and angle follows a uniform law.
r = abs(randn(N,1)*source.lambda);

% initial direction
theta = 2*pi*rand(N,1);
d = theta;

% initial positions in cartesian coordinates
x = r.*cos(theta);
y = r.*sin(theta);

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
if physics.acoustics && isfield(material,'v')
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
P = struct('N', N, ...              % number of particles
           'r', r, ...              % cylindrical coordinates
           'theta', theta, ...              
           'x', x, ...              % cartesian coordinates
           'y', y, ...              
           'd', d, ...              % propagation direction
           'p', p, ...              % polarization (used only in elasticity)
           'meanFreePath', mfp, ... % mean free path
           'v', v, ...              % propagation velocity
           't', t );                % current time for the particle

