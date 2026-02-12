function obs = radiativeTransfer( geometry, source, material, observation )

% physics
d = geometry.dimension;
acoustics = material.acoustics;
material.timeSteps = 0;
if ~isfield(geometry,'frame')
    geometry.frame = 'spherical';
end

% discretization in packets of particles  (for optimal vectorization)
Npk = 1e5; % size of packets 
Np = ceil(source.numberParticles/Npk); % number of packets

% initialize observation structure
[ obs, E, bins, ibins, vals, Nt, t , d1, d2 ] = ...
         initializeObservation( geometry, acoustics, observation, Np*Npk );

% prepare scattering cross sections 
material = prepareSigma( material, d );      

times = zeros(1, Np); 
tic; % Start overall timer

% loop on packages of particles
for ip = 1:Np
    tic; % Start iteration timer

    % initialize particles
    P = initializeParticle( Npk, d, acoustics, source );

    % loop on time
    for it = 2:Nt

        % propagate particles
        if material.timeSteps==0
            P = propagateParticleSmallDt( material, geometry, P, t(it) );
        elseif material.timeSteps==1
            P = propagateParticle( material, P, t(it) );
        end

        % observe energies (as a function of [Psi x t])
        E(:,:,it,:) = E(:,:,it,:) + observeTime( geometry, acoustics, ...
                                      P.x, P.p, P.dir, bins, ibins, vals );

    % end of loop on time
    end

    % end of loop on packages

    % Store iteration time
    times(ip) = toc;
    
    % Calculate estimates every 10 iterations or at start
    if ip == 1 || mod(ip, 10) == 0
        avg_time = mean(times(1:ip)); 
        elapsed = toc; 
        remaining_iterations = Np - ip;
        estimated_remaining = avg_time * remaining_iterations;
        estimated_total = elapsed + estimated_remaining;
        % Display progress
        fprintf('Iteration %d/%d (%.1f%%)\n', ip, Np, ip/Np*100);
        fprintf('  Elapsed: %.2f s\n', elapsed);
        fprintf('  Estimated remaining: %.2f s (%.2f min)\n', ...
                estimated_remaining, estimated_remaining/60);
        fprintf('  Estimated total: %.2f s (%.2f min)\n\n', ...
                estimated_total, estimated_total/60);
    end
end

fprintf('Total elapsed time: %.2f s (%.2f min)\n\n', ...
                sum(times), sum(times)/60);

% energy density as a function of [x1 x2 t]
obs.energyDensity = (1./(d1'*(d2*obs.N))).*double(E);

% energy as a function of [t]
obs.energy = squeeze(tensorprod(d1,reshape(tensorprod(d2,obs.energyDensity,2,2),[length(d1) Nt 1+~acoustics]),2,1));
