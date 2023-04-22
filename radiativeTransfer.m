function obs = radiativeTransfer( source, material, observation, geometry )

% construct set of virtual sources
Lmax = material.v*max(observation.time);
[ns,posS,Rmax] = virtualSources( geometry, source.position, Lmax);
observation.r = linspace(0,Rmax,ceil(Rmax/observation.dr));
 
% compute solution in full space
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% construct superposition of total energy for different sources
dx = ceil( geometry.size(1)/observation.dr );
obs.boxX = linspace( 0, geometry.size(1), dx );
dz = ceil( geometry.size(3)/observation.dr );
obs.boxZ = linspace( -geometry.size(3), 0, dz );
[boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);
E = zeros(numel(boxx),obs.Nt);
for i1 = 1:ns
    r = sqrt((boxx-posS(i1,1)).^2+(posS(1,2)-posS(i1,2)).^2+(boxz-posS(i1,3)).^2);
    E = E + interp1(obs.r',obs.energyDensity,r(:),'linear',0);
end
E = permute(reshape(E,length(obs.boxZ),length(obs.boxX),obs.Nt),[2 1 3]);
obs.energyDensityBox = E;

% construct directional energy at sensors
obs.sensors = observation.sensors;
Nss = size(obs.sensors,1);
obs.psi2pi = [0 obs.psi 2*pi-obs.psi(end:-1:1) 2*pi];
obs.energyDirectional = zeros(2*(obs.Npsi+1),obs.Nt,Nss);
for i1 = 1:Nss
    for i2 = 1:ns
        X = obs.sensors(i1,:) - posS(i2,:);
        r = sqrt(sum(X.^2));
        rxz = sqrt(sum(X([1 3]).^2));
        Es = squeeze( interp1( obs.r', permute( obs.energy,[2 1 3]), r ));
        Es = [ Es(1,:); Es; Es(end,:); Es(end:-1:1,:) ];
        theta = acos(X(1)/rxz);
        if X(3)<0; theta = 2*pi-theta; end
        psi = mod( obs.psi2pi - theta, 2*pi );
        obs.energyDirectional(:,:,i1) = obs.energyDirectional(:,:,i1) ...
                                         + interp1( obs.psi2pi', Es, psi );
    end
end

end
