function obs = radiativeTransferUnbounded( geometry, source, material, observation )

% physics
d = geometry.dimension;
acoustics = material.acoustics;

% discretization in packets of particles  (for optimal vectorization)
Npk = 1e5; % size of packets 
Np = ceil(source.numberParticles/Npk); % number of packets

% initialize observation structure
[ obs, E, bins, ibins, vals, Nt, t , d1, d2 ] = ...
         initializeObservation( geometry, acoustics, observation, Np*Npk );

% prepare scattering cross sections 
material = prepareSigma( material, d );      

% loop on packages of particles
for ip = 1:Np

    % initialize particles
    P = initializeParticle( Npk, d, acoustics, source );

    % loop on time
    for it = 2:Nt

        % propagate particles
        if material.timeSteps==0
            P = propagateParticleSmallDt( material, geometry.bnd, P, t(it) );
        elseif material.timeSteps==1
            P = propagateParticle( material, P, t(it) );
        end

        % observe energies (as a function of [Psi x t])
        E(:,:,it,:) = E(:,:,it,:) + observeTime( geometry, acoustics, ...
                                      P.x, P.p, P.dir, bins, ibins, vals );

    % end of loop on time
    end

% end of loop on packages
end

% energy density as a function of [x1 x2 t]
obs.energyDensity = (1./(d1'*(d2*obs.N))).*double(E);

% energy as a function of [t]
obs.energy = tensorprod(d1,squeeze(tensorprod(d2,obs.energyDensity,2,2)),2,1);
