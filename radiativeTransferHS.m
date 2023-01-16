function obs = radiativeTransferHS( source, material, observation, geometry )

% compute solution in full space
observation.movie = true;
obs = radiativeTransfer( source, material, observation );

% compute solution in lower half-space (1)
z0 = geometry.sourcePosition;
if z0>0; error('source should be in lower half-space'); end
obs.x = obs.grid;
obs.z = obs.grid+z0;
E1 = obs.gridEnergy;

% compute solution in upper half-space (2)
z2 = obs.grid-z0;
E2 = permute(obs.gridEnergy,[2 1 3]);
E2 = interp1(z2,E2,obs.z,'linear',0);
E2 = permute(E2,[2 1 3]);
if strcmpi(geometry.boundaryCondition,'dirichlet')
    E2 = -E2;
end

% merge both solutions
obs.zmin = min(z2);
obs.energyHS = E1+E2;