function obs = radiativeTransfer( physics, source, material, observation, geometry )

% construct set of virtual sources
Lmax = material.v*max(observation.time);
[ns,posS,signS,Rmax] = virtualSources( geometry, source.position, Lmax, physics.dimension);
observation.sensors = linspace(0,Rmax,ceil(Rmax/observation.dx));
 
% compute solution in full space
obs = radiativeTransferUnbounded( physics, source, material, observation );

% construct superposition of solutions for different sources
dx = ceil( geometry.size(1)/observation.dx );
obs.boxX = linspace( 0, geometry.size(1), dx );
dz = ceil( geometry.size(3)/observation.dx );
obs.boxZ = linspace( -geometry.size(3), 0, dz );
[boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);
E = zeros(numel(boxx),obs.Nt);
for i1 = 1:ns
    r = sqrt((boxx-posS(i1,1)).^2+(boxz-posS(i1,3)).^2);
    E = E + signS(i1)*interp1(obs.x',obs.energyDensity,r(:),'linear',0);
end
E = permute(reshape(E,length(obs.boxZ),length(obs.boxX),obs.Nt),[2 1 3]);
obs.energyDensityBox = E;
end

% - idea to be able to consider non-isotropic initial conditions: perform the 
% simulation with initial direction uniform, and keep in memory the initial
% direction of each particle. And, when performing the histogram, keep one
% more dimension in the histograms that identicates exactly that initial
% direction. Then, when considering a non-isotropic initial condition, it
% is only necessary to perform a multiplication along that dimension with
% the desired distribution for the initial direction
%
% - create a routine that would take as input the PSD directly of
% mechanical/acoustical parameters and return the scattering cross sections
% in the acoustical and elasticity cases (Shahram me les a déjà envoyées ?)
%
% - create routines that allow to choose loadings and forces depending
% directly on the corresponding quantities for the wave equation
%
