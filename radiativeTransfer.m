function obs = radiativeTransfer( geometry, source, material, observation )

simulationTimer = tic;

% physics
d = geometry.dimension;
acoustics = material.acoustics;

if ~isfield(geometry,'frame')
    geometry.frame = 'spherical'; % default frame is spherical
end

% Check that at most one frequency is provided.
if ~isempty(material.Frequency) && ~isscalar(material.Frequency)
    error('material.Frequency must be a scalar.');
end

% Check the inputs needed for deterministic acoustic attenuation.
useAcousticAttenuation = acoustics && ~isempty(material.Q) && ...
    isfinite(material.Q(1));
if useAcousticAttenuation
    if material.Q(1) <= 0
        error('The acoustic quality factor Q must be positive.');
    end
    if isempty(material.Frequency)
        error(['A finite quality factor Q was provided, but material.Frequency is empty. ', ...
               'Please define material.Frequency before using Q-based attenuation.']);
    end
end

% Check if user provided a single Q value for dissipation of both P and S waves
if isprop(material,'Q') && ~isempty(material.Q) && ~material.acoustics
    if isscalar(material.Q) && isfinite(material.Q)
        warning(['Only one quality factor was provided for an elastic simulation. ', ...
                 'The same value Q = %g will be used for both P and S waves.'], material.Q);
    end
end

% Checks if Parallel Toolbox exists
hasParallelToolbox = ~isempty(ver('parallel'));

% discretization in packets of particles for optimal vectorization
Npk = 1e5;
Ntotal = source.numberParticles;
Np = ceil(Ntotal/Npk); % number of packets
particlesPerPacket = Npk * ones(1,Np);
particlesPerPacket(end) = Ntotal - Npk*(Np-1);

% initialize observation structure
[ obs, E, bins, ibins, vals, Nt, t , d1, d2 ] = ...
    initializeObservation( geometry, acoustics, observation, Ntotal );

% prepare scattering cross sections
material = MaterialClass.prepareSigma( material, d );

printSimulationSummary(geometry,source,material,Ntotal,Np, ...
    hasParallelToolbox);

% Pre-calculate properties of E so workers know what to allocate
szE = size(E);
classE = class(E);

dt = mean(diff(t));
t_start = t(1);

if hasParallelToolbox
    % =====================================================================
    % PARALLEL EXECUTION (PARFOR)
    % =====================================================================
    % This sends 'material' and 'geometry' to workers ONLY ONCE.
    cMaterial = parallel.pool.Constant(material);
    cGeometry = parallel.pool.Constant(geometry);

    % Progress tracking via DataQueue
    parNp        = Np;
    parDone      = 0;
    parTstart    = tic;
    parQueue     = parallel.pool.DataQueue;
    afterEach(parQueue, @updateParProgress_);

    % loop on packages of particles
    parfor ip = 1:Np

        Npacket = particlesPerPacket(ip);

        % Retrieve local copies from the Constant wrapper
        matLocal = cMaterial.Value;
        geoLocal = cGeometry.Value;

        % Local accumulator
        E_local = zeros(szE, classE);

        % initialize particles
        P = initializeParticle( Npacket, d, acoustics, source );

        % The propagation loop starts at it = 2. Record the normalized
        % initial condition in the first frame (t = 0+, after injection).
        E_local(:,:,1,:) = observeTime( geoLocal, acoustics, ...
            P.x, P.p, P.dir, bins, ibins, vals, P.alive );

        % loop on time
        if matLocal.timeSteps == 0
            for it = 2:Nt
                t_val = t_start + (it-1)*dt;
                P = propagateParticleSmallDt( matLocal, geoLocal, P, t_val );
                E_local(:,:,it,:) = E_local(:,:,it,:) + observeTime( geoLocal, acoustics, ...
                    P.x, P.p, P.dir, bins, ibins, vals, P.alive );
            end
        else
            for it = 2:Nt
                t_val = t_start + (it-1)*dt;
                P = propagateParticle( matLocal, P, t_val );
                E_local(:,:,it,:) = E_local(:,:,it,:) + observeTime( geoLocal, acoustics, ...
                    P.x, P.p, P.dir, bins, ibins, vals, P.alive );
            end
        end

        % Add the worker's local result to the main variable
        E = E + E_local;

        % Signal completion of this packet to the progress tracker
        send(parQueue, true);

    end

else
    % =====================================================================
    % SERIAL EXECUTION (FOR)
    % =====================================================================
    times = zeros(1, Np);
    serialTimer = tic;

    % loop on packages of particles
    for ip = 1:Np
        packetTimer = tic;

        Npacket = particlesPerPacket(ip);

        % initialize particles
        P = initializeParticle( Npacket, d, acoustics, source );

        % The propagation loop starts at it = 2. Record the normalized
        % initial condition in the first frame (t = 0+, after injection).
        E(:,:,1,:) = E(:,:,1,:) + observeTime( geometry, acoustics, ...
            P.x, P.p, P.dir, bins, ibins, vals, P.alive );

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
                P.x, P.p, P.dir, bins, ibins, vals, P.alive );

            % end of loop on time
        end

        % end of loop on packages

        % Store iteration time
        times(ip) = toc(packetTimer);

        % Calculate estimates every 10 iterations or at start
        if ip == 1 || mod(ip, 10) == 0
            avg_time = mean(times(1:ip));
            elapsed = toc(serialTimer);
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

        % end of loop on packages
    end
end

% Apply intrinsic attenuation deterministically to acoustic energies.
E = double(E);
if useAcousticAttenuation
    omega = 2*pi*material.Frequency;
    attenuation = reshape(exp(-omega*t/material.Q(1)),1,1,Nt,1);
    E = E .* attenuation;
end

% energy density as a function of [x1 x2 t]
obs.energyDensity = (1./(d1'*(d2*obs.N))).*E;

% correc the normalizations to get an energy density
if obs.d == 3
    if isequal(ibins, [1 2])
        obs.energyDensity = obs.energyDensity / obs.dz;
    elseif isequal(ibins, [1 3])
        obs.energyDensity = obs.energyDensity / obs.dy;
    elseif isequal(ibins, [1 4])
        obs.energyDensity = obs.energyDensity .* obs.dpsi / (obs.dy*obs.dz);
    elseif isequal(ibins, [2 4])
        obs.energyDensity = obs.energyDensity .* obs.dpsi / (obs.dx*obs.dz);
    elseif isequal(ibins, [3 4])
        obs.energyDensity = obs.energyDensity .* obs.dpsi / (obs.dx*obs.dy);
    elseif isequal(ibins, [2 3])
        obs.energyDensity = obs.energyDensity / obs.dx;
    else
        error('Unsupported bin combination in 3D');
    end
else
    if isequal(ibins, [1 4])
        obs.energyDensity = obs.energyDensity .* obs.dpsi / obs.dy;
    elseif isequal(ibins, [2 4])
        obs.energyDensity = obs.energyDensity .* obs.dpsi / obs.dx;
    elseif ~isequal(ibins, [1 2])
        error('Unsupported bin combination in 2D');
    end
end

% energy as a function of [t]
obs.energy = squeeze(sum(sum(E,1),2)) / obs.N;

fprintf('Radiative Transfer Monte Carlo simulation completed in %.2f s.\n\n', ...
    toc(simulationTimer));

    % -----------------------------------------------------------------
    % Nested function: called by afterEach on the DataQueue.
    % Shares parNp, parDone, parTstart with the parent workspace.
    function updateParProgress_( ~ )
        parDone  = parDone + 1;
        elapsed  = toc(parTstart);
        printEvery = max(1, floor(parNp / 20)); % ~5 % steps
        if parDone == 1 || mod(parDone, printEvery) == 0 || parDone == parNp
            eta = elapsed / parDone * (parNp - parDone);
            fprintf('  Packet %d/%d (%.0f%%) | Elapsed: %.1fs | ETA: %.1fs (%.1fmin)\n', ...
                parDone, parNp, parDone/parNp*100, elapsed, eta, eta/60);
        end
    end

end % radiativeTransfer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printSimulationSummary(geometry,source,material,Ntotal,Npackets, ...
        hasParallelToolbox)

physics = 'elastic';
if material.acoustics
    physics = 'acoustic';
end

scatteringModel = material.scatteringModel;
if strcmpi(scatteringModel,'parameterFluctuations')
    scatteringModel = 'material-parameter fluctuations';
end

correlationModel = material.SpectralLaw;
if strcmpi(material.scatteringModel,'polycrystal') && ~isempty(material.TPCF)
    correlationModel = material.TPCF.model;
elseif isempty(correlationModel) && ~isempty(material.sigma)
    correlationModel = 'user-defined DSCS';
end
if isempty(correlationModel)
    correlationModel = 'not specified';
end

sourceType = 'point';
if isfield(source,'type') && ~isempty(source.type)
    sourceType = char(source.type);
end
if ~material.acoustics
    polarization = 'P';
    if isfield(source,'polarization') && ~isempty(source.polarization)
        polarization = char(source.polarization);
    end
    sourceType = sprintf('%s, initially %s polarized',sourceType,polarization);
end

sourcePosition = [0 0 0];
if isfield(source,'position') && ~isempty(source.position)
    sourcePosition = source.position;
end

sourceDirection = 'isotropic';
if isfield(source,'direction') && ~isempty(source.direction)
    if isnumeric(source.direction)
        sourceDirection = mat2str(source.direction,6);
    else
        sourceDirection = char(source.direction);
        if strcmpi(sourceDirection,'upper')
            sourceDirection = 'upper hemisphere';
        end
    end
end

fprintf('\n------------------------------------------------------------\n');
fprintf('Radiative Transfer Monte Carlo simulation\n');
fprintf('  Problem               : %d-D %s, %s coordinates\n', ...
    geometry.dimension,physics,geometry.frame);
fprintf('  Scattering model      : %s\n',scatteringModel);
fprintf('  Correlation structure : %s\n',material.correlationStructure);
fprintf('  Correlation model     : %s\n',correlationModel);
if isempty(material.Frequency)
    fprintf('  Frequency             : not specified\n');
else
    fprintf('  Frequency             : %.6g Hz\n',material.Frequency);
end
if material.acoustics
    fprintf('  Wave velocity         : %.6g m/s\n',material.v);
else
    fprintf('  Wave velocities       : Vp = %.6g m/s, Vs = %.6g m/s\n', ...
        material.vp,material.vs);
end
fprintf('  Source type           : %s\n',sourceType);
fprintf('  Source position       : %s m\n',mat2str(sourcePosition,6));
fprintf('  Source direction      : %s\n',sourceDirection);
if isfield(source,'lambda') && ~isempty(source.lambda)
    fprintf('  Source spatial width  : %.6g m\n',source.lambda);
end
fprintf('  Particles             : %d in %d packets\n',Ntotal,Npackets);

if ~isfield(geometry,'bnd') || isempty(geometry.bnd)
    fprintf('  Boundaries            : none (unbounded medium)\n');
else
    boundaryType = 'reflective';
    if isfield(geometry.bnd,'type')
        specifiedTypes = {geometry.bnd.type};
        specifiedTypes = specifiedTypes(~cellfun('isempty',specifiedTypes));
        if ~isempty(specifiedTypes) && all(strcmpi(specifiedTypes{1},specifiedTypes))
            boundaryType = lower(specifiedTypes{1});
        elseif ~isempty(specifiedTypes)
            boundaryType = 'mixed';
        end
    end
    fprintf('  Boundaries            : %d %s\n', ...
        numel(geometry.bnd),boundaryType);
end

if isempty(material.Q) || all(isinf(material.Q(:)))
    fprintf('  Intrinsic attenuation : none\n');
elseif material.acoustics
    fprintf('  Intrinsic attenuation : Q = %.6g\n',material.Q(1));
elseif isscalar(material.Q)
    fprintf('  Intrinsic attenuation : Qp = Qs = %.6g\n',material.Q);
else
    fprintf('  Intrinsic attenuation : Qp = %.6g, Qs = %.6g\n', ...
        material.Q(1),material.Q(2));
end

if all(material.Sigma(:) == 0)
    fprintf('No volume scattering detected: the medium is homogeneous.\n');
    fprintf(['Particles will propagate ballistically; boundary interactions ', ...
        'and intrinsic attenuation remain active.\n']);
end

if material.acoustics
    fprintf('  Acoustic mean free time          : %.6g s\n', ...
        material.meanFreeTime);
    fprintf('  Acoustic mean free path          : %.6g m\n', ...
        material.meanFreePath);
    fprintf('  Acoustic transport mean free time: %.6g s\n', ...
        material.transportMeanFreeTime);
    fprintf('  Acoustic transport mean free path: %.6g m\n', ...
        material.transportMeanFreePath);
else
    fprintf('  P and S mean free times          : %.6g s, %.6g s\n', ...
        material.meanFreeTime);
    fprintf('  P and S mean free paths          : %.6g m, %.6g m\n', ...
        material.meanFreePath);
    fprintf('  P and S transport mean free times: %.6g s, %.6g s\n', ...
        material.transportMeanFreeTime);
    fprintf('  P and S transport mean free paths: %.6g m, %.6g m\n', ...
        material.transportMeanFreePath);
end
fprintf('------------------------------------------------------------\n\n');

end
