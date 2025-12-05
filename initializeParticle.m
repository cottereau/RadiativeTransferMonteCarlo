function P = initializeParticle( N, d, acoustics, source )

% default source
if ~isfield(source,'type')
    source.type = 'point';
end
if ~isfield(source,'position')
    source.position = [0 0 0];
end
if d==2 && length(source.position)==2
    source.position = [source.position 0];
end

% initial position of each particle
switch source.type

    % point source
    case 'point'
        % radius follows a Gaussian law with standard deviation lambda, and
        % angle follows a uniform law. If radial option is used, a 
        % user-defined function (of the radius) is used instead.
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
        % random rotation of the perp vector in 3D
        if d==3
            alpha = 2*pi*rand(N,1);
            perp = perp.*cos(alpha) + cross(dir,perp,2).*sin(alpha);
        end

        % initial positions in cartesian coordinates
        x = source.position + ...
            r.*[cos(theta0).*sin(phi0) sin(theta0).*sin(phi0) cos(phi0)];

    % plane waves
    case 'plane'
        extent = zeros(N,3);
        extent(:,setdiff(1:3,source.direction)) = repmat(source.extent,[N 1]);
        extent(:,source.direction) = randn(N,1)*source.lambda/2;
        x = source.position + (rand(N,3)-0.5).*extent;
        dir = zeros(N,3);
        dir(:,abs(source.direction)) = sign(source.direction);
        perp = zeros(N,3);
        theta = rand(N,1)*2*pi;
        perp(:,setdiff(1:3,source.direction)) = [cos(theta) sin(theta)];

end

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

% [~,ind] = unique(sort(cdf(ind)));
%invcdf = griddedInterpolant(sort(cdf(ind)),sort(xth(ind)),'linear','nearest');

% Sort and remove duplicates
[unique_cdf, idx] = unique(sort(cdf(ind)));
ds = sort(xth(ind));
unique_xth = ds(idx);

% Ensure monotonicity
if ~issorted(unique_cdf)
    [unique_cdf, sort_idx] = sort(unique_cdf);
    unique_xth = unique_xth(sort_idx);
end

% Create the interpolant
invcdf = griddedInterpolant(unique_cdf, unique_xth, 'linear', 'nearest');

end
