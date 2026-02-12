function obs = radiativeTransferUnbounded( geometry, source, material, observation )

% physics
d = geometry.dimension;
acoustics = material.acoustics;

if ~isfield(geometry,'frame')
    geometry.frame = 'spherical'; % default frame is spherical
end

% Checks if Parallel Toolbox exists
hasParallelToolbox = ~isempty(ver('parallel'));

% discretization in packets of particles for optimal vectorization
Npk = 1e5;
Np = ceil(source.numberParticles/Npk); % number of packets

% initialize observation structure
[ obs, E, bins, ibins, vals, Nt, t , d1, d2 ] = ...
         initializeObservation( geometry, acoustics, observation, Np*Npk );

% prepare scattering cross sections 
material = prepareSigma( material, d );      

% Pre-calculate properties of E so workers know what to allocate
szE = size(E);
classE = class(E);

dt = mean(diff(t));
t_start = t(1);

% Start the global timer
tStart = tic;

if hasParallelToolbox
    % =====================================================================
    % PARALLEL EXECUTION (PARFOR)
    % =====================================================================
    fprintf('Parallel Computing Toolbox detected. Running on multiple workers\n');
    
    % This sends 'material' and 'geometry' to workers ONLY ONCE.
    cMaterial = parallel.pool.Constant(material);
    cGeometry = parallel.pool.Constant(geometry);

    % loop on packages of particles
    parfor ip = 1:Np

        % Retrieve local copies from the Constant wrapper
        matLocal = cMaterial.Value;
        geoLocal = cGeometry.Value;

        % Local accumulator
        E_local = zeros(szE, classE);
    
        % initialize particles
        P = initializeParticle( Npk, d, acoustics, source );
        
        % loop on time
        if matLocal.timeSteps == 0
            for it = 2:Nt
                t_val = t_start + (it-1)*dt;
                P = propagateParticleSmallDt( matLocal, geoLocal, P, t_val );
                E_local(:,:,it,:) = E_local(:,:,it,:) + observeTime( geoLocal, acoustics, ...
                    P.x, P.p, P.dir, bins, ibins, vals );
            end
        else
            for it = 2:Nt
                t_val = t_start + (it-1)*dt;
                P = propagateParticle( matLocal, P, t_val );
                E_local(:,:,it,:) = E_local(:,:,it,:) + observeTime( geoLocal, acoustics, ...
                    P.x, P.p, P.dir, bins, ibins, vals );
            end
        end
        
        % Add the worker's local result to the main variable
        E = E + E_local;
        
    end

else
    % =====================================================================
    % SERIAL EXECUTION (FOR)
    % =====================================================================
    fprintf('Parallel Computing Toolbox NOT found. Running in serial mode\n');

    % loop on packages of particles
    for ip = 1:Np
    
        % initialize particles
        P = initializeParticle( Npk, d, acoustics, source );
        
        % NOTE: In serial, we do NOT need E_local. 
        % We can write directly to E, saving memory and overhead.
    
        % loop on time
        for it = 2:Nt
    
            % Calculate time locally
            t_val = t_start + (it-1)*dt; 
    
            % propagate particles
            if material.timeSteps==0
                P = propagateParticleSmallDt( material, geometry, P, t_val );
            elseif material.timeSteps==1
                P = propagateParticle( material, P, t_val );
            end
    
            % observe energies - DIRECT UPDATE
            E(:,:,it,:) = E(:,:,it,:) + observeTime( geometry, acoustics, ...
                                       P.x, P.p, P.dir, bins, ibins, vals );
    
        % end of loop on time
        end
        
        % Progress reporting is reliable in serial loops
        if mod(ip, 10) == 0 || ip == Np
             tElapsed = toc(tStart);
             fprintf('Package %d / %d finished. (Elapsed: %.1fs)\n', ip, Np, tElapsed);
        end
        
    % end of loop on packages
    end

end

% energy density as a function of [x1 x2 t]
obs.energyDensity = (1./(d1'*(d2*obs.N))).*double(E);

% energy as a function of [t]
obs.energy = squeeze(tensorprod(d1,reshape(tensorprod(d2,obs.energyDensity,2,2),[length(d1) Nt 1+~acoustics]),2,1));
