function P = initializeParticle( N, d, acoustics, source )

% default source
if ~isfield(source,'position')
    source.position = [0 0 0];
end
if d==2 && length(source.position)==2
    source.position = [source.position 0];
end

% initial position of each particle: radius follows a Gaussian law with
% standard deviation lambda, and angle follows a uniform law.
% if radial option is used, a user-defined function (of the radius) is used
% instead
if ~isfield(source,'radial') || isempty(source.radial)
    r = abs(randn(N,1)*source.lambda/2);
else
    Rmax = source.radial.GridVectors{1}(end);
    invcdfsource = inverseCDF(source.radial,d,Rmax);
    r = invcdfsource(rand(N,1));
end

% initial angles and directions
if d==2
    phi0 = (pi/2)*ones(N,1);
    phi = (pi/2)*ones(N,1);
elseif d==3
    phi0 = acos(-1+2*rand(N,1));
    phi = acos(-1+2*rand(N,1));
end
theta0 = 2*pi*rand(N,1);
theta = 2*pi*rand(N,1);
if isfield(source,'direction') && strcmp(source.direction,'outgoing')
    theta = theta0;
    phi = phi0;
elseif isfield(source,'direction') && strcmp(source.direction,'plane')
    theta = zeros(N,1);
    phi = (pi/2)*ones(N,1);
end
dir = [cos(theta).*sin(phi) sin(theta).*sin(phi) cos(phi)];
perp = [-sin(theta) cos(theta) zeros(N,1)];

% initial positions in cartesian coordinates
x = source.position + ...
              r.*[cos(theta0).*sin(phi0) sin(theta0).*sin(phi0) cos(phi0)];

% current time for the particle
t = zeros(N,1);

% initial polarisation of each particle
% in acoustics, this variable is unused (and set always to true)
% in elastics, true corresponds to P waves (default), and false to S waves
p = true(N,1);
if ~acoustics && isfield(source,'polarization') && source.polarization=='S'
    p = false(N,1);
end

% initialize structure
P = struct( 'd', d, ...                 % dimension of the problem
            'N', N, ...                 % number of particles
            'x', x, ...                 % cartesian coordinates
            'dir', dir, ...             % direction of propagation
            'perp', perp, ...           % orthogonal to direction of propagation
            'p', p, ...                 % polarization (used only in elasticity)
            't', t );                   % current time for the particle

end

% compute the cumulative distribution radial function corresponding to a
% given (positive) function to draw randomly from it
function invcdf = inverseCDF(f,d,Rmax)
Nth = 10000;
xth = linspace(0,Rmax,Nth);
intF = trapz(xth,(xth.^(d-1)).*f(xth));
pdf = (xth.^(d-1)).*f(xth)/intF;
cdf = cumsum(pdf)*mean(diff(xth));
cdf(1) = 0;
cdf(end) = 1;
ind = find(diff(cdf)>0);
ind = unique([ind ind+1]);
invcdf = griddedInterpolant(cdf(ind),xth(ind),'linear','nearest');
end
