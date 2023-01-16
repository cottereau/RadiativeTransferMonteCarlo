function P = initializeParticle(N,source,material)

% initialization
P = struct('t', zeros(N,1), ...
           'N', N );

% initial position of each particle: radius follows a Gaussian law with
% standard deviation lambda, and angle follows a uniform law.
P.r = abs(randn(N,1)*source.lambda);
%P.theta = zeros(N,1);
%P.x = P.r*ones(N,1);
%P.y = zeros(N,1);


% initial polarisation of each particle
% in acoustics, this variable is unused (and set always to true)
% in elastics, true corresponds to P waves, and false to S waves
P.p = true(N,1);
if isfield(source,'polarization') && source.polarization=='S'
    P.p = false(N,1);
end

% initial direction of propagation of each particle (angle)
% P.theta = zeros(N,1);
% P.costheta = cos(P.theta);
P.costheta = ones(N,1);

% material velocity of the background for each particle
if isfield(material,'v')
    P.v = material.v*ones(N,1);
elseif isfield(material,'vp') && isfield(material,'vs')
    P.v = material.vs*ones(N,1);
    P.v(p) = material.vp;
else
    error('unknown velocity of the background')
end
% meanFreePath for each particle
if isfield(material,'meanFreePath')
    P.meanFreePath = material.meanFreePath*ones(N,1);
elseif isfield(material,'meanFreePathP') && isfield(material,'meanFreePathS')
    P.meanFreePath = material.meanFreePathS*ones(N,1);
    P.meanFreePath(p) = material.meanFreePathP;
else
    error('unknown meanFreePath')
end
