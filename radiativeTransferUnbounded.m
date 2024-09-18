function obs = radiativeTransferUnbounded( d, source, material, observation )

% physics
acoustics = material.acoustics;

% discretization in packets of particles  (for optimal vectorization)
Npk = 1e4; % size of packets 
Np = ceil(source.numberParticles/Npk); % number of packets

% initialize observation structure
[ obs, E, bins, ibins, vals, Nt, t , d1, d2 ] = ...
                initializeObservation( d, acoustics, observation, Np*Npk );

% prepare scattering cross sections 
material = prepareSigma( material, d );      

% loop on packages of particles
parfor ip = 1:Np
    % initialize particles
    P = initializeParticle( Npk, d, acoustics, source );
    x = zeros( Npk, 3, Nt ); x(:,:,1) = P.x;
    p = true( Npk, Nt ); p(:,1) = P.p;
    dir = zeros( Npk, 3, Nt ); dir(:,:,1) = P.dir;

    % loop on time
    for it = 2:Nt

        % propagate particles
        P = propagateParticle( material, P, t(it) );

        % observe energies (as a function of [Psi x t])
        x(:,:,it) = P.x;
        p(:,it) = P.p;
        dir(:,:,it) = P.dir;

    % end of loop on time
    end
    
    % aggregate results
    E = E + observeTime( d, acoustics, x, p, dir, bins, ibins, vals );

% end of loop on packages
end
obs.energyDensity = (1./(d1'*(d2*obs.N))).*double(E);

% energy as a function of [t]
obs.energy = tensorprod(d1,squeeze(tensorprod(d2,obs.energyDensity,2,2)),2,1);
