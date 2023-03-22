function obs = radiativeTransfer( physics, source, material, observation, geometry )

% construct set of virtual sources
Lmax = material.v*max(observation.time);
[ns,posS,signS,Rmax] = virtualSources( geometry, source.position, Lmax);
observation.sensors = linspace(0,Rmax,ceil(Rmax/observation.dx));

% compute solution in full space
obs = radiativeTransferUnbounded( physics, source, material, observation );

% construct superposition of solutions for different sources
[~,obs.boxX,obs.boxZ] = computBox(geometry,source.position,observation.dx);
[boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);
E = zeros(numel(boxx),obs.Nt);
for i1 = 1:ns
    r = sqrt((boxx-posS(i1,1)).^2+(boxz-posS(i1,3)).^2);
    E = E + signS(i1)*interp1(obs.x',obs.energyDensity,r(:),'linear',0);
end
E = permute(reshape(E,length(obs.boxZ),length(obs.boxX),obs.Nt),[2 1 3]);
obs.energyDensityBox = E;
end

