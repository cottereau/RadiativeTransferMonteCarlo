function obs = radiativeTransfer( source, material, observation, geometry )

% construct set of virtual sources
Lmax = material.v*max(observation.time);
[ ns, posSources, Rmax ] = virtualSources( geometry, source.position, Lmax);
observation.r = linspace(0,Rmax,ceil(Rmax/observation.dr));
 
% compute solution in full space
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
obs.nSources = ns;
obs.positionSources = posSources;

% preparing plotting the energy
dx = ceil( geometry.size(1)/observation.dr );
obs.boxX = linspace( 0, geometry.size(1), dx );
dz = ceil( geometry.size(3)/observation.dr );
obs.boxZ = linspace( -geometry.size(3), 0, dz );
[boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);

% local amplifications for slabs and boxes
maskin = ones(numel(boxx),1);
maskout = ones(numel(boxx),1);
if strcmp(geometry.type,'box') || strcmp(geometry.type,'slab')
    zin = [-geometry.size(3) 0]+[-1 1]*source.position(3);
    maskin = maskin.*(1+exp(-(boxz(:)-zin(1)).^2/source.lambda^2));
    maskin = maskin.*(1+exp(-(boxz(:)-zin(2)).^2/source.lambda^2));
    zout = [-geometry.size(3) 0];
    maskout = maskout.*(1+exp(-(boxz(:)-zout(1)).^2/source.lambda^2));
    maskout = maskout.*(1+exp(-(boxz(:)-zout(2)).^2/source.lambda^2));
end
if strcmp(geometry.type,'box')
    xin = [geometry.size(1) 0]+[-1 1]*source.position(1);
    maskin = maskin.*(1+exp(-(boxx(:)-xin(1)).^2/source.lambda^2));
    maskin = maskin.*(1+exp(-(boxx(:)-xin(2)).^2/source.lambda^2));
    xout = [geometry.size(1) 0];
    maskout = maskout.*(1+exp(-(boxx(:)-xout(1)).^2/source.lambda^2));
    maskout = maskout.*(1+exp(-(boxx(:)-xout(2)).^2/source.lambda^2));
end
energyCoherent = squeeze(tensorprod(obs.dpsi',obs.energy(:,:,:,1),1));
energyIncoherent = squeeze(tensorprod(obs.dpsi',obs.energy(:,:,:,2),1));

% construct superposition of total energy for different sources
E = zeros(numel(boxx),obs.Nt);
for i1 = 1:ns
    x = boxx-posSources(i1,1);
    y = posSources(1,2)-posSources(i1,2);
    z = boxz-posSources(i1,3);
    r = sqrt(x.^2+y.^2+z.^2);
    Ecoherent = interp1(obs.r',energyCoherent,r(:),'linear',0);
    Ecoherent = Ecoherent .* repmat(maskin.*maskout,[1 obs.Nt]);
    Eincoherent = interp1(obs.r',energyIncoherent,r(:),'linear',0);
    Eincoherent = Eincoherent .* repmat(maskout,[1 obs.Nt]);
    E = E + Ecoherent + Eincoherent;
end
E = permute(reshape(E,length(obs.boxZ),length(obs.boxX),obs.Nt),[2 1 3]);
obs.energyDensityBox = E;

end
