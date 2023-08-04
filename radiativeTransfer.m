function obs = radiativeTransfer( source, material, observation, geometry )

% constants
lambda = source.lambda;
posS0 = source.position;
type = geometry.type;
siz = geometry.size;
d = geometry.dimension;
v = material.v;

% construct set of virtual sources
Lmax = material.v*max(observation.time);
[ ns, posS, Rmax ] = virtualSources( geometry, source.position, Lmax);
observation.r = linspace(0,Rmax,ceil(Rmax/observation.dr));
 
% compute solution in full space
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
energyIncoherent = squeeze(sum(obs.energy(:,:,:,2),2));
obs.nSources = ns;
obs.positionSources = posS;

% preparing plotting the energy
dx = ceil( siz(1)/observation.dr );
obs.boxX = linspace( 0, siz(1), dx );
dz = ceil( siz(3)/observation.dr );
obs.boxZ = linspace( -siz(3), 0, dz );
[boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);

% local amplifications at boundaries for slabs and boxes
amp = ones(size(boxx));
if strcmp(type,'halfspace') || strcmp(type,'slab') || strcmp(type,'box')
    amp = amp .* (1 + exp(-(boxz/lambda).^2));
end
if  strcmp(type,'slab') || strcmp(type,'box')
    amp = amp .* (1 + exp(-((boxz+siz(3))/lambda).^2));
end
if  strcmp(type,'box')
    amp = amp .* (1 + exp(-(boxx/lambda).^2));
    amp = amp .* (1 + exp(-((boxx-siz(1))/lambda).^2));
    amp = amp .* (1 + exp(-(posS0(2)/lambda).^2));
    amp = amp .* (1 + exp(-((posS0(2)-siz(3))/lambda).^2));
end
amp = amp(:);

% Superposition of incoherent energy coming from each source
E = zeros(numel(boxx),obs.Nt);
for i1 = 1:ns
    x = boxx(:)-posS(i1,1);
    y = posS(1,2)-posS(i1,2);
    z = boxz(:)-posS(i1,3);
    r = sqrt(x.^2+y.^2+z.^2);
    Eincoherent = interp1(obs.r',energyIncoherent,r(:),'linear',0);
    E = E + Eincoherent.*amp;
end

% add coherent energy
Ecoherent = coherentInABox(obs.Ec,boxx(:),posS(1,2),boxz(:),posS,obs.t,d,lambda,v);
E = E + Ecoherent;

% reformatting
E = permute(reshape(E,length(obs.boxZ),length(obs.boxX),obs.Nt),[2 1 3]);
obs.energyDensityBox = E;

end
