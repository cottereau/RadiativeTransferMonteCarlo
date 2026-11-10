function obs = radiativeTransfer( geometry, source, material, observation, opts )
%RADIATIVETRANSFER  Monte Carlo solution of the radiative transfer equation.
%   obs = radiativeTransfer(geometry, source, material, observation)
%   obs = radiativeTransfer(geometry, source, material, observation, opts)
%
%   Monochromatic unit-energy response to a source event. material is a
%   MaterialClass object, or a cell array {mat1, mat2, ...} for several regions
%   separated by geometry.bnd(k).mat = [i j] interfaces (particles start in
%   material{1}, which also supplies acoustics, Frequency and Q).
%   source is a source structure (point, plane, disk, or a physical source from
%   sourcePointForce / sourceMomentTensor / sourceAcoustic) or a cell array of
%   sources with .fraction fields.
%
%   opts (all optional):
%     eventTime            source emission time; earlier observation times are empty (0)
%     packetSize           particles per packet (1e5)
%     scatteringAngleNodes nodes of the scattering-angle CDFs (1e6)
%     polarized            unbounded 3-D Born pure-state transport (false)
%     maxPolarizedParticles particle budget of the polarized path (1e6)
%     binAverage           record at each output time the average over its time bin
%                          (midpoints between output times, clipped at the event);
%                          each particle is observed at its own uniform random time
%                          in the bin. Requires polarized transport or timeSteps = 1 (false)
%
%   obs.energy is the captured energy [time, mode] (modes P, S for elastic waves);
%   obs.aliveEnergy also counts the surviving energy outside the receiver grid.

simulationTimer = tic;
% multi-material: material = {mat1, mat2, ...} with geometry.bnd(k).mat = [i j]
% interfaces (propagateParticleMultiMaterial); particles start in material{1};
% material{1} also supplies acoustics/Frequency/Q/timeSteps below.
matArray = {};
if iscell(material)
    matArray = material;
    material = matArray{1};
end

% options
if nargin < 5, opts = struct; end
if ~isstruct(opts) || ~isscalar(opts), error('RTE:Options','Options must be a scalar structure.'); end
defaults = struct('eventTime',0,'packetSize',1e5,'scatteringAngleNodes',1e6, ...
    'polarized',false,'maxPolarizedParticles',1e6,'binAverage',false);
names = fieldnames(defaults);
if any(~ismember(fieldnames(opts),names)), error('RTE:Options','Unknown solver option.'); end
for j = 1:numel(names)
    if ~isfield(opts,names{j}), opts.(names{j}) = defaults.(names{j}); end
end
validateattributes(opts.eventTime,{'numeric'},{'real','finite','scalar'});
validateattributes(opts.packetSize,{'numeric'},{'real','finite','scalar','integer','positive','<=',flintmax});
validateattributes(opts.scatteringAngleNodes,{'numeric'},{'real','finite','scalar','integer','>=',16});
if ~islogical(opts.binAverage) || ~isscalar(opts.binAverage)
    error('RTE:Options','binAverage must be a logical scalar.');
end
if ~islogical(opts.polarized) || ~isscalar(opts.polarized)
    error('Fonte:PolarizedOptions','polarized must be a logical scalar.');
end
validateattributes(opts.maxPolarizedParticles,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'});

% physics
d = geometry.dimension;
acoustics = material.acoustics;
validateattributes(d,{'numeric'},{'scalar','integer','>=',2,'<=',3});
if ~isfield(geometry,'frame')
    geometry.frame = 'spherical'; % default frame is spherical
end
hasBoundaries = isfield(geometry,'bnd') && ~isempty(geometry.bnd);
if opts.polarized && (d~=3 || acoustics || material.d~=3 || ~isempty(matArray) || ...
        ~strcmp(material.scatteringModel,'parameterFluctuations') || ...
        ~strcmp(material.correlationStructure,'isotropic') || ...
        ~isa(material.Phi,'function_handle') || material.timeSteps~=1 || hasBoundaries || ...
        (~isempty(material.Q) && any(isfinite(material.Q))))
    error('Fonte:PolarizedScope', ['Polarized transport currently requires an unbounded ', ...
        '3D elastic single material with timeSteps=1 and no intrinsic attenuation.']);
end

% sources
if iscell(source), sourceList = source; else, sourceList = {source}; end
if isempty(sourceList), error('Fonte:InvalidSource','At least one source is required.'); end
for j = 1:numel(sourceList)
    s = sourceList{j};
    if ~isstruct(s) || ~isscalar(s), error('Fonte:InvalidSource','Sources must be scalar structures.'); end
    if isfield(s,'type') && any(strcmp(s.type,{'point_force','moment_tensor','acoustic'}))
        physicalSourceForContext(s,d,acoustics,material);
    end
end
if isfield(sourceList{1},'numberParticles')
    N = sourceList{1}.numberParticles;
elseif iscell(source)
    N = 1e5;
else
    error('Fonte:MissingSourceField','Source requires numberParticles.');
end
validateattributes(N,{'numeric'},{'real','finite','scalar','integer','positive','<=',flintmax});
N = double(N);
if opts.polarized && N>opts.maxPolarizedParticles
    error('Fonte:PolarizedBudget','Polarized particle count exceeds the explicit budget.');
end

% material checks
if ~isempty(material.Frequency)
    validateattributes(material.Frequency,{'numeric'},{'real','finite','scalar','nonnegative'});
end
validateattributes(material.timeSteps,{'numeric'},{'scalar','integer','>=',0,'<=',1});
if ~isempty(material.Q)
    validateattributes(material.Q,{'numeric'},{'real','positive','vector','nonnan'});
    if ~(isscalar(material.Q) || (~acoustics && numel(material.Q)==2))
        error('RTE:QualityFactor','Q must be scalar or an elastic [QP QS] pair.');
    end
end
useAcousticAttenuation = acoustics && ~isempty(material.Q) && isfinite(material.Q(1));
if ~isempty(material.Q) && any(isfinite(material.Q)) && isempty(material.Frequency)
    error('RTE:FrequencyRequired','Finite Q requires a material frequency.');
end
if ~acoustics && isscalar(material.Q) && isfinite(material.Q)
    warning(['Only one quality factor was provided for an elastic simulation. ', ...
             'The same value Q = %g will be used for both P and S waves.'], material.Q);
end
if isempty(matArray) && material.timeSteps==1
    if hasBoundaries
        error('RTE:UnsupportedBoundary','Boundary propagation requires timeSteps=0.');
    end
    if ~acoustics && ~isempty(material.Q) && any(isfinite(material.Q))
        error('RTE:UnsupportedAttenuation','Elastic attenuation requires timeSteps=0.');
    end
end
if opts.binAverage && (~isempty(matArray) || (~opts.polarized && material.timeSteps~=1))
    error('RTE:BinAverage','binAverage requires polarized transport or timeSteps=1 (single material).');
end

% source populations allocated once, then cut into packets
plan = sourcePacketPlan(source,N,d,acoustics);
packets = zeros(0,2);                                  % [plan group, particle count]
for g = 1:numel(plan)
    remaining = plan(g).count;
    while remaining > 0
        count = min(double(opts.packetSize),remaining);
        packets(end+1,:) = [g count]; %#ok<AGROW>
        remaining = remaining - count;
    end
end
Np = size(packets,1);

% observation
[obs,E,bins,ibins,vals,Nt,t] = initializeObservation(geometry,acoustics,observation,N);
E = double(E);  nmodes = 1+~acoustics;
age = t - double(opts.eventTime);
aliveCounts = zeros(Nt,nmodes);  launched = zeros(1,nmodes);
% time bins for binAverage: midpoints between output times, symmetric at the ends
a = age(:);  edges = [a(1); (a(1:end-1)+a(2:end))/2; a(end)];
if Nt>1, edges([1 end]) = [a(1)-(a(2)-a(1))/2; a(end)+(a(end)-a(end-1))/2]; end

% prepare scattering cross sections (only needed after the event)
prepared = any(age>0);
if prepared
    if isempty(matArray)
        material = MaterialClass.prepareSigma(material,d,opts.scatteringAngleNodes);
    else
        for k = 1:numel(matArray)
            matArray{k} = prepareRegion(matArray{k},d,opts.scatteringAngleNodes);
        end
        material = matArray{1};
    end
end

hasParallelToolbox = ~isempty(ver('parallel'));
printSimulationSummary(geometry,source,material,N,Np,hasParallelToolbox);

% no source particle is visible before the common emission event
if any(age>=0)
    szE = size(E);
    if hasParallelToolbox
        % =================================================================
        % PARALLEL EXECUTION (PARFOR)
        % =================================================================
        % Workers get the client path, so they resolve the same classes even
        % when another folder of the saved path defines a MaterialClass.
        pool = gcp;
        wait(parfevalOnAll(pool, @path, 0, path));

        % 'material' and 'geometry' are sent to the workers only once.
        cMaterial = parallel.pool.Constant(material);
        cGeometry = parallel.pool.Constant(geometry);

        % progress tracking via DataQueue
        parNp     = Np;
        parDone   = 0;
        parTstart = tic;
        parQueue  = parallel.pool.DataQueue;
        afterEach(parQueue, @updateParProgress_);

        parfor ip = 1:Np
            [Ek,Ak,Lk] = runPacket(plan(packets(ip,1)).source, packets(ip,2), d, acoustics, ...
                cMaterial.Value, matArray, cGeometry.Value, bins, ibins, vals, age, edges, opts, szE);
            E = E + Ek;
            aliveCounts = aliveCounts + Ak;
            launched = launched + Lk;
            send(parQueue, true);
        end
    else
        % =================================================================
        % SERIAL EXECUTION (FOR)
        % =================================================================
        times = zeros(1, Np);
        serialTimer = tic;
        for ip = 1:Np
            packetTimer = tic;
            [Ek,Ak,Lk] = runPacket(plan(packets(ip,1)).source, packets(ip,2), d, acoustics, ...
                material, matArray, geometry, bins, ibins, vals, age, edges, opts, szE);
            E = E + Ek;
            aliveCounts = aliveCounts + Ak;
            launched = launched + Lk;

            % progress estimate every 10 packets
            times(ip) = toc(packetTimer);
            if ip == 1 || mod(ip, 10) == 0
                elapsed = toc(serialTimer);
                estimated_remaining = mean(times(1:ip)) * (Np - ip);
                fprintf('Iteration %d/%d (%.1f%%)\n', ip, Np, ip/Np*100);
                fprintf('  Elapsed: %.2f s\n', elapsed);
                fprintf('  Estimated remaining: %.2f s (%.2f min)\n', ...
                    estimated_remaining, estimated_remaining/60);
                fprintf('  Estimated total: %.2f s (%.2f min)\n\n', ...
                    elapsed + estimated_remaining, (elapsed + estimated_remaining)/60);
            end
        end
    end
end

% intrinsic attenuation applied deterministically to acoustic energies
if useAcousticAttenuation
    attenuation = exp(-2*pi*material.Frequency*max(age,0)/material.Q(1));
    E = E.*reshape(attenuation,1,1,Nt,1);
    aliveCounts = aliveCounts.*attenuation(:);
end

% Direction-integrated cells use the spatial volume; resolved directions also
% divide by their angle/cosine width. Infinite spatial intervals are marginals.
obs.energyDensity = E./(N*obs.cellMeasure);
obs.energy = reshape(sum(sum(E,1),2),Nt,nmodes)/N;
obs.capturedEnergy = obs.energy;
obs.aliveEnergy = aliveCounts/N;
obs.outsideObservationEnergy = max(0,sum(obs.aliveEnergy,2)-sum(obs.energy,2));
obs.emittedEnergy = double(age(:)>=0);
obs.eventTime = double(opts.eventTime);
obs.age = age;
obs.launchModeFraction = launched/N;
obs.unitInjectedEnergy = true;
obs.scatteringPrepared = prepared;
obs.options = opts;
obs.modeLabels = {'P','S'};
if acoustics, obs.modeLabels = {'A'}; end

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

% -------------------------------------------------------------------------
function [E,aliveCounts,launched] = runPacket(src, count, d, acoustics, material, matArray, ...
        geometry, bins, ibins, vals, age, edges, opts, szE)
% One packet: draw the particles, advance them through the output times and
% record what each time sees.
Nt = numel(age);  nmodes = 1+~acoustics;
E = zeros(szE);  aliveCounts = zeros(Nt,nmodes);  launched = zeros(1,nmodes);
P = initializeParticle(count,d,acoustics,src);
launched(1) = nnz(P.p);
if ~acoustics, launched(2) = nnz(~P.p); end
for it = 1:Nt
    if age(it) < 0, continue; end
    if opts.binAverage
        lo = max(edges(it),0);
        P = advance(material,matArray,geometry,P,lo,opts);
        Q = advance(material,matArray,geometry,P,lo+rand(count,1)*(edges(it+1)-lo),opts);
    else
        P = advance(material,matArray,geometry,P,age(it),opts);  Q = P;
    end
    E(:,:,it,:) = E(:,:,it,:) + observeTime(geometry,acoustics, ...
        Q.x,Q.p,Q.dir,bins,ibins,vals,Q.alive);
    aliveCounts(it,1) = nnz(Q.alive & Q.p);
    if ~acoustics, aliveCounts(it,2) = nnz(Q.alive & ~Q.p); end
end
end

% -------------------------------------------------------------------------
function P = advance(material,matArray,geometry,P,T,opts)
% Propagate to time T (scalar, or per particle for timeSteps=1 / polarized).
if all(T<=0), return; end
if ~isempty(matArray)
    P = propagateParticleMultiMaterial(matArray,geometry,P,T);
elseif opts.polarized
    P = propagatePolarizedParticle(material,P,T);
elseif material.timeSteps==0
    P = propagateParticleSmallDt(material,geometry,P,T);
else
    P = propagateParticle(material,P,T);
end
end

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

if isempty(material.Sigma)
    fprintf('  Scattering            : not prepared (no output time after the event)\n');
elseif all(material.Sigma(:) == 0)
    fprintf('No volume scattering detected: the medium is homogeneous.\n');
    fprintf(['Particles will propagate ballistically; boundary interactions ', ...
        'and intrinsic attenuation remain active.\n']);
end

if isempty(material.meanFreeTime)
    % nothing prepared
elseif material.acoustics
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

% -------------------------------------------------------------------------
function mat = prepareRegion( mat, d, Nth )
% A region with zero coefficients of variation (e.g. a homogeneous water
% layer) does not scatter: prepareSigma cannot build a CDF from sigma == 0,
% so it gets an infinite mean free time directly.
if isempty(mat.sigma) && mat.acoustics && ~any(mat.coefficients_of_variation)
    mat.vp = mat.v; mat.vs = 0;
    mat.Sigma = 0; mat.Sigmapr = 0;
    mat.meanFreeTime = Inf; mat.meanFreePath = Inf;
    mat.invcdf = @(u) zeros(size(u)); % never drawn
else
    mat = MaterialClass.prepareSigma( mat, d, Nth );
end
end
