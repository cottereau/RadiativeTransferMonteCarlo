function out = radiativeTransferBroadband(geometry,spectrum,materialFactory,observation,opts)
%RADIATIVETRANSFERBROADBAND Physical P/S energy from monochromatic responses.
% spectrum = sourceEnergySpectrum(source,wave); factory(f_Hz) must return a
% fresh MaterialClass for every active band. For an acoustic source
% (sourceAcoustic) the factory returns acoustic materials with v = c, and
% the single mode is reported as {'A'}. All frequencies are in Hz at
% the factory boundary, rad/s in bands. No spectral renormalization occurs.
% opts: bands, bandOptions, solverOptions, seed, batches, maxRuns,
% maxTotalParticles, keepBandResponses, responseFrequencies. Defaults cap
% work to tiny runs. responseFrequencies (Hz) computes the unit responses
% only at these nodes and interpolates them piecewise-linearly in omega
% (clamped) onto the source quadrature. Since the source weights enter
% linearly, this equals transport at the nodes with the band energies
% lumped by the hat functions: E_k = sum_j L_k(omega_j) w_j E'(omega_j).
% The lumped energies are out.responseBands; the source spectrum is still
% integrated on the fine grid, the response is smooth in omega.
if nargin<5, opts=struct; end
if ~isstruct(opts) || ~isscalar(opts), error('Fonte:BroadbandOptions','Options must be a scalar structure.'); end
defaults=struct('bands',[],'bandOptions',struct,'solverOptions',struct, ...
    'seed',1704,'batches',1,'maxRuns',8,'maxTotalParticles',100, ...
    'keepBandResponses',false,'responseFrequencies',[]);
names=fieldnames(defaults);
if any(~ismember(fieldnames(opts),names)), error('Fonte:BroadbandOptions','Unknown broadband option.'); end
for k=1:numel(names)
    if ~isfield(opts,names{k}), opts.(names{k})=defaults.(names{k}); end
end
validateattributes(opts.seed,{'numeric'},{'real','finite','scalar','integer','nonnegative','<=',2^32-1});
validateattributes(opts.batches,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(opts.maxRuns,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(opts.maxTotalParticles,{'numeric'},{'real','finite','scalar','integer','positive'});
if ~isscalar(opts.keepBandResponses) || ~islogical(opts.keepBandResponses)
    error('Fonte:BroadbandOptions','keepBandResponses must be logical.');
end
if ~isa(materialFactory,'function_handle') || ~isstruct(spectrum) || ...
        ~isfield(spectrum,'source') || ~isfield(spectrum,'P') || ...
        ~isfield(spectrum,'S') || ~isfield(spectrum,'eventTime')
    error('Fonte:BroadbandInput','Use a physical source spectrum and material factory.');
end
if isempty(opts.bands)
    bands=sourceFrequencyBands(spectrum,opts.bandOptions);
else
    bands=opts.bands;
end
required={'omega','frequencyHz','weightOmega','energyP','energyS','eventTime'};
if ~isstruct(bands) || ~isscalar(bands) || ~all(isfield(bands,required))
    error('Fonte:InvalidBands','Bands need nodes, weights, modal energies and event time.');
end
omega=double(bands.omega(:)); frequency=double(bands.frequencyHz(:));
weight=double(bands.weightOmega(:)); p=double(bands.energyP(:)); s=double(bands.energyS(:));
n=numel(omega);
if n<1 || any([numel(frequency) numel(weight) numel(p) numel(s)]~=n) || ...
        any(~isfinite([omega;frequency;weight;p;s])) || ...
        any(omega<=0) || any(weight<=0) || any(p<0) || any(s<0) || ...
        any(diff(omega)<=0) || any(abs(frequency-omega/(2*pi))>1e-12*frequency) || ...
        bands.eventTime~=spectrum.eventTime
    error('Fonte:InvalidBands','Invalid positive-frequency quadrature or event time.');
end
expectedP=weight.*spectrum.P(omega); expectedS=weight.*spectrum.S(omega);
if any(abs(p-expectedP)>1e-9*max(expectedP,realmin)) || ...
        any(abs(s-expectedS)>1e-9*max(expectedS,realmin))
    error('Fonte:InvalidBands','Band energies do not match this source spectrum.');
end
responseBands=bands;
if ~isempty(opts.responseFrequencies)
    validateattributes(opts.responseFrequencies,{'numeric'},{'real','finite','positive','vector'});
    node=2*pi*double(opts.responseFrequencies(:));
    if any(diff(node)<=0), error('Fonte:BroadbandOptions','responseFrequencies must increase strictly.'); end
    L=hatWeights(omega,node);
    p=L'*p; s=L'*s; omega=node; frequency=node/(2*pi); n=numel(node);
    responseBands=struct('omega',omega,'frequencyHz',frequency,'energyP',p,'energyS',s, ...
        'eventTime',bands.eventTime,'interpolation','piecewise linear in omega, clamped');
end
src=spectrum.source;
if ~isfield(src,'numberParticles') || ~isfield(src,'type') || ...
        ~any(strcmp(src.type,{'point_force','moment_tensor','acoustic'}))
    error('Fonte:BroadbandInput','Spectrum needs a physical source with numberParticles.');
end
acoustic=strcmp(src.type,'acoustic');
validateattributes(src.numberParticles,{'numeric'},{'real','finite','scalar','integer','positive'});
if ~isstruct(opts.solverOptions) || ~isscalar(opts.solverOptions) || ...
        isfield(opts.solverOptions,'eventTime')
    error('Fonte:BroadbandOptions','solverOptions must be a struct; eventTime comes from the spectrum.');
end
solverOptions=opts.solverOptions; solverOptions.eventTime=bands.eventTime;
active=[p>0 s>0]; activeBands=any(active,2);
runCount=nnz(active)*opts.batches;
if runCount==0, error('Fonte:ZeroBandEnergy','No positive-energy channels.'); end
if runCount>opts.maxRuns || runCount*double(src.numberParticles)>opts.maxTotalParticles
    error('Fonte:BroadbandBudget','Requested %d runs and %g particles exceed the explicit work budget.', ...
        runCount,runCount*double(src.numberParticles));
end
oldRng=rng; cleanup=onCleanup(@()rng(oldRng)); %#ok<NASGU>
meta=repmat(struct('frequencyHz',[],'omega',[],'materialClass',[], ...
    'Sigma',[],'meanFreeTime',[],'seeds',[]),n,1);
if opts.keepBandResponses, kept=cell(n,2); end
priorMaterial=[]; out=[]; varDensity=[]; varEnergy=[];
for j=1:n
    if ~activeBands(j), continue; end
    material=materialFactory(frequency(j));
    if ~isa(material,'MaterialClass') || ~isscalar(material) || ...
            isempty(material.Frequency) || ~isscalar(material.Frequency) || ...
            abs(double(material.Frequency)-frequency(j))>1e-12*frequency(j) || ...
            (~isempty(priorMaterial) && material==priorMaterial)
        error('Fonte:MaterialFactory','Factory must return a fresh material at the requested frequency.');
    end
    if acoustic
        if ~isequal(material.acoustics,true) || ...
                ~isnumeric(material.v) || ~isscalar(material.v) || ...
                abs(material.v-src.c)>1e-12*src.c || ...
                (~isempty(material.rho) && abs(material.rho-src.rho)>1e-12*src.rho)
            error('Fonte:MaterialMismatch','Each material must be acoustic and match the source c and rho.');
        end
    elseif ~isequal(material.acoustics,false) || ...
            ~isnumeric(material.vp) || ~isscalar(material.vp) || ...
            ~isnumeric(material.vs) || ~isscalar(material.vs) || ...
            ~isnumeric(material.rho) || ~isscalar(material.rho) || ...
            any(abs([material.vp material.vs material.rho]-[src.vp src.vs src.rho]) > ...
            1e-12*[src.vp src.vs src.rho])
        error('Fonte:MaterialMismatch','Each material must match the source background vp, vs and rho.');
    end
    priorMaterial=material;
    meta(j).frequencyHz=frequency(j); meta(j).omega=omega(j);
    meta(j).materialClass=class(material);
    bandResponses=cell(1,2); seeds=nan(2,opts.batches);
    channelEnergy=[p(j) s(j)];
    for beta=1:2
        if ~active(j,beta), continue; end
        conditional=src; conditional.launchMode='P';
        if beta==2, conditional.launchMode='S'; end
        meanDensity=[]; meanEnergy=[]; meanAlive=[];
        for b=1:opts.batches
            seed=mod(double(opts.seed)+double((j-1)*2*opts.batches+(beta-1)*opts.batches+b-1),2^32);
            seeds(beta,b)=seed; rng(seed,'twister');
            r=radiativeTransfer(geometry,conditional,material,observation,solverOptions);
            if b==1
                meanDensity=zeros(size(r.energyDensity)); meanEnergy=zeros(size(r.energy));
                meanAlive=zeros(size(r.aliveEnergy));
                m2Density=meanDensity; m2Energy=meanEnergy;
            end
            d=r.energyDensity-meanDensity; meanDensity=meanDensity+d/b;
            m2Density=m2Density+d.*(r.energyDensity-meanDensity);
            d=r.energy-meanEnergy; meanEnergy=meanEnergy+d/b;
            m2Energy=m2Energy+d.*(r.energy-meanEnergy);
            meanAlive=meanAlive+(r.aliveEnergy-meanAlive)/b;
        end
        r.energyDensity=meanDensity; r.energy=meanEnergy; r.aliveEnergy=meanAlive;
        bandResponses{beta}=r;
        if opts.batches>1
            vd=channelEnergy(beta)^2*m2Density/(opts.batches*(opts.batches-1));
            ve=channelEnergy(beta)^2*m2Energy/(opts.batches*(opts.batches-1));
            if isempty(varDensity), varDensity=zeros(size(vd)); varEnergy=zeros(size(ve)); end
            varDensity=varDensity+vd; varEnergy=varEnergy+ve;
        end
        if opts.keepBandResponses, kept{j,beta}=r; end
    end
    meta(j).Sigma=material.Sigma; meta(j).meanFreeTime=material.meanFreeTime;
    meta(j).seeds=seeds;
    oneBand=struct('energyP',p(j),'energyS',s(j),'eventTime',bands.eventTime);
    part=composeSourceBands(oneBand,bandResponses);
    if isempty(out)
        out=part;
    else
        if ~isequal(part.t,out.t) || ~isequal(part.cellMeasure,out.cellMeasure) || ...
                ~isequal(size(part.energyDensity),size(out.energyDensity))
            error('Fonte:ResponseMismatch','Band responses have different observation grids.');
        end
        out.energyDensity=out.energyDensity+part.energyDensity;
        out.energy=out.energy+part.energy;
        out.aliveEnergy=out.aliveEnergy+part.aliveEnergy;
    end
end
out.capturedEnergy=out.energy;
out.outsideObservationEnergy=max(0,sum(out.aliveEnergy,2)-sum(out.energy,2));
out.bandEnergy=[sum(p) sum(s)];
out.emittedEnergy=sum(out.bandEnergy)*double(out.t(:)>=bands.eventTime);
out.source=spectrum.source; out.waveform=spectrum.waveform;
out.bands=bands; out.responseBands=responseBands; out.materialMetadata=meta;
out.geometry=geometry; out.observation=observation;
out.options=opts; out.solverOptions=solverOptions;
out.runCount=runCount; out.particleCount=runCount*double(src.numberParticles);
if acoustic
    out.transportModel='acoustic scalar Born scattering';
elseif isfield(solverOptions,'polarized') && solverOptions.polarized
    out.transportModel='3D polarized Born scattering, exact rates, pure-state unraveling';
else
    out.transportModel='legacy scalar-scattering approximation';
end
if opts.batches>1
    out.mcStandardError=struct('energyDensity',sqrt(varDensity),'energy',sqrt(varEnergy));
else
    out.mcStandardError=struct('energyDensity',[],'energy',[]);
end
if opts.keepBandResponses, out.bandResponses=kept; end
end

function L=hatWeights(x,node)
% L(j,k): weight of response node k at fine frequency x(j); rows sum to one.
L=zeros(numel(x),numel(node));
if isscalar(node), L(:)=1; return; end
x=min(max(x,node(1)),node(end));
k=min(discretize(x,node),numel(node)-1);
t=(x-node(k))./(node(k+1)-node(k)); j=(1:numel(x))';
L(sub2ind(size(L),j,k))=1-t;
L(sub2ind(size(L),j,k+1))=t;
end
