function P = initializeParticle(N,source,material)

% dimension
d = material.dimension;

% initial position of each particle: radius follows a Gaussian law with
% standard deviation lambda, and angle follows a uniform law.
r = abs(randn(N,1)*source.lambda);

% initial angles and directions
if d==2
    theta = (pi/2)*ones(N,1);
elseif d==3
    theta = 2*pi*rand(N,1);
end
phi = 2*pi*rand(N,1);
dth = theta;
dphi = phi;

% initial positions in cartesian coordinates
x = r.*sin(theta).*cos(phi);
y = r.*cos(theta).*sin(phi);
z = r.*cos(theta);

% initial polarisation of each particle
% in acoustics, this variable is unused (and set always to true)
% in elastics, true corresponds to P waves, and false to S waves
p = true(N,1);
if isfield(source,'polarization') && source.polarization=='S'
    p = false(N,1);
end

% current time for the particle
t = zeros(N,1);

% material velocity of the background for each particle
if isfield(material,'v')
    v = material.v*ones(N,1);
elseif isfield(material,'vp') && isfield(material,'vs')
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
           'r', r, ...              % spherical coordinates
           'theta', theta, ...      %  theta in [0  pi]        
           'phi', phi, ...          %  phi   in [0 2pi]
           'x', x, ...              % cartesian coordinates
           'y', y, ...              
           'z', z, ...              
           'dth', dth, ...          % propagation direction
           'dphi', dphi, ...        
           'p', p, ...              % polarization (used only in elasticity)
           'meanFreePath', mfp, ... % mean free path
           'v', v, ...              % propagation velocity
           't', t );                % current time for the particle

