function obs = radiativeTransfer( source, material, observation, geometry )

% construct set of virtual sources
Lmax = material.v*max(observation.time);
[ ns, posSources, Rmax ] = virtualSources( geometry, source.position, Lmax);
observation.r = linspace(0,Rmax,ceil(Rmax/observation.dr));
 
% compute solution in full space
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
obs.nSources = ns;
obs.positionSources = posSources;

% construct superposition of total energy for different sources
dx = ceil( geometry.size(1)/observation.dr );
obs.boxX = linspace( 0, geometry.size(1), dx );
dz = ceil( geometry.size(3)/observation.dr );
obs.boxZ = linspace( -geometry.size(3), 0, dz );
[boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);
E = zeros(numel(boxx),obs.Nt);
for i1 = 1:ns
    x = boxx-posSources(i1,1);
    y = posSources(1,2)-posSources(i1,2);
    z = boxz-posSources(i1,3);
    r = sqrt(x.^2+y.^2+z.^2);
    E = E + interp1(obs.r',obs.energyDensity,r(:),'linear',0);
end
E = permute(reshape(E,length(obs.boxZ),length(obs.boxX),obs.Nt),[2 1 3]);
obs.energyDensityBox = E;


end
