function obs = radiativeTransfer( source, material, observation, geometry )

% construct set of virtual sources
Lmax = material.v*max(observation.time);
[ns,posS,signS,Rmax] = virtualSources( geometry, source.position, Lmax);
observation.sensors = linspace(0,Rmax,ceil(Rmax/observation.dx));
 
% compute solution in full space
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% construct superposition of solutions for different sources
dx = ceil( geometry.size(1)/observation.dx );
obs.boxX = linspace( 0, geometry.size(1), dx );
dz = ceil( geometry.size(3)/observation.dx );
obs.boxZ = linspace( -geometry.size(3), 0, dz );
[boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);
E = zeros(numel(boxx),obs.Nt);
for i1 = 1:ns
    r = sqrt((boxx-posS(i1,1)).^2+(posS(1,2)-posS(i1,2)).^2+(boxz-posS(i1,3)).^2);
    E = E + signS(i1)*interp1(obs.x',obs.energyDensity,r(:),'linear',0);
end
E = permute(reshape(E,length(obs.boxZ),length(obs.boxX),obs.Nt),[2 1 3]);
obs.energyDensityBox = E;

% construct directional energy at sensors
obs.sensors = observation.sensors;
i1 = 1; i2 = 1;
r = sum((obs.sensors(i1,:)-posS(i2,:)).^2);
end
