function comparison = compareRickerBallisticSource(mechanism,opts)
%COMPARERICKERBALLISTICSOURCE Light RTE versus analytic ballistic cell energy.
% The analytic reference is the event-source, zero-scattering limit of the
% Ricker or supplied sampled source: E_beta times its angular pattern on
% the shell r=c_beta*(t-t0). This is not the finite-width coherent pulse.
if nargin<2, opts=struct; end
if ~isstruct(opts) || ~isscalar(opts), error('Fonte:ExampleOptions','Options must be a scalar structure.'); end
defaults=struct('numberParticles',8,'seed',1801,'vp',4000,'vs',2000, ...
    'rho',2500,'amplitude',1,'tensor',[],'waveform',[], ...
    'sampledTime',[],'sampledValues',[]);
names=fieldnames(defaults);
if any(~ismember(fieldnames(opts),names)), error('Fonte:ExampleOptions','Unknown comparison option.'); end
for j=1:numel(names)
    if ~isfield(opts,names{j}), opts.(names{j})=defaults.(names{j}); end
end
validateattributes(opts.numberParticles,{'numeric'},{'real','finite','scalar','integer','>=',1,'<=',8});
validateattributes(opts.seed,{'numeric'},{'real','finite','scalar','integer','nonnegative','<=',2^32-1});
oldPath=path; oldRng=rng;
cleanup=onCleanup(@()restore(oldPath,oldRng)); %#ok<NASGU>
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root,'-begin'); AddRTMCLib(root,false); rng(opts.seed,'twister');
caseOptions=struct('vp',opts.vp,'vs',opts.vs,'rho',opts.rho, ...
    'amplitude',opts.amplitude,'numberParticles',opts.numberParticles, ...
    'tensor',opts.tensor);
example=fonteBenchmarkCase(mechanism,'R',caseOptions);
source=example.source;
if ~isempty(opts.sampledTime) || ~isempty(opts.sampledValues)
    if ~isempty(opts.waveform)
        error('Fonte:ExampleOptions','Choose either a waveform or sampled time/values.');
    end
    opts.waveform=sourceTimeFunction('sampled', ...
        struct('time',opts.sampledTime,'values',opts.sampledValues));
end
if isempty(opts.waveform)
    wave=example.waveform; spec=example.spectrum; caseID=example.caseID;
else
    wave=opts.waveform; spec=sourceEnergySpectrum(source,wave);
    caseID=[upper(char(mechanism)) '-SAMPLED'];
end
geometry=struct('dimension',3,'frame','spherical');
observation=struct('x',[1000 2500 5000 10000 20000 40000], 'y',[-pi pi], ...
    'z',[-pi/2 -pi/6 pi/6 pi/2], 'directions',[0 pi], ...
    'time',spec.eventTime+[1.01 2.53 5.11]);
solverOptions=struct('eventTime',spec.eventTime,'packetSize',opts.numberParticles, ...
    'scatteringAngleNodes',16,'polarized',true, ...
    'maxPolarizedParticles',8);
fractions=angularFractions(source,observation.z);
energies=spec.energy;
caseResults=cell(1,2); labels={'P','S'};
for beta=1:2
    if energies(beta)==0, continue; end
    conditioned=source; conditioned.launchMode=labels{beta};
    material=zeroMaterial(opts.vp,opts.vs,opts.rho,10);
    response=radiativeTransfer(geometry,conditioned,material, ...
        observation,solverOptions);
    if any(abs(response.energy(:,3-beta))>1e-12)
        error('Fonte:ExampleMode','Zero scattering unexpectedly converted modes.');
    end
    measured=response.energyDensity(:,:,:,beta).*response.cellMeasure.*energies(beta);
    expected=zeros(size(measured));
    finitePulse=zeros(size(measured));
    speed=opts.vp; if beta==2, speed=opts.vs; end
    for it=1:numel(response.t)
        radius=speed*(response.t(it)-spec.eventTime);
        ir=find(radius>=observation.x(1:end-1) & radius<observation.x(2:end),1);
        if isempty(ir), error('Fonte:ExampleWindow','Ballistic front left the receiver grid.'); end
        expected(ir,:,it)=energies(beta)*fractions(beta,:);
        elsewhere=measured(:,:,it); elsewhere(ir,:)=0;
        if any(abs(elsewhere(:))>1e-11*energies(beta))
            error('Fonte:ExampleFront','RTE energy left the analytic ballistic shell.');
        end
        for jr=1:numel(observation.x)-1
            lo=response.t(it)-observation.x(jr+1)/speed;
            hi=response.t(it)-observation.x(jr)/speed;
            pulseFraction=pulseEnergyFraction(wave,spec.derivativeOrder,lo,hi);
            finitePulse(jr,:,it)=energies(beta)*pulseFraction*fractions(beta,:);
        end
    end
    captured=sum(sum(measured,1),2);
    if any(abs(captured(:)-energies(beta))>1e-11*energies(beta))
        error('Fonte:ExampleEnergy','RTE did not preserve all-space modal energy.');
    end
    countError=energies(beta)*sqrt(fractions(beta,:).*(1-fractions(beta,:))/opts.numberParticles);
    relativeRMSE=sqrt(mean(((measured(:)-expected(:))/energies(beta)).^2));
    caseResults{beta}=struct('mode',labels{beta},'energyJ',energies(beta), ...
        'measuredCellEnergyJ',measured,'analyticCellEnergyJ',expected, ...
        'finitePulseFarFieldCellEnergyJ',finitePulse, ...
        'measuredDensity',measured./response.cellMeasure, ...
        'analyticDensity',expected./response.cellMeasure, ...
        'finitePulseFarFieldDensity',finitePulse./response.cellMeasure, ...
        'angularFraction',fractions(beta,:), ...
        'binomialStandardErrorJ',countError, ...
        'relativeCellRMSE',relativeRMSE,'unitResponse',response);
end
comparison=struct('caseID',caseID,'source',source,'waveform',wave, ...
    'spectrum',spec,'geometry',geometry,'observation',observation, ...
    'P',caseResults{1},'S',caseResults{2},'expectedRatioSP', ...
    energies(2)/energies(1),'numberParticlesPerMode',opts.numberParticles, ...
    'analyticModel','full-spectrum event-source ballistic cell energy', ...
    'note','RTE is the event-source model. The separate finite-pulse curve uses far-field radiated power and an interpolated sampled waveform.');
end

function fraction=pulseEnergyFraction(wave,m,a,b)
% Far-field shell energy fraction: integral_{a}^{b} |s^(m)(tau)|^2 dtau.
if strcmp(wave.kind,'sampled')
    t=double(wave.parameters.time(:));
    a=max(a,t(1)); b=min(b,t(end));
    if b<=a, fraction=0; return; end
    edges=unique([a;t(t>a & t<b);b]);
else
    f=wave.parameters.f0; center=wave.parameters.delay;
    a=max(a,center-8/(pi*f)); b=min(b,center+8/(pi*f));
    if b<=a, fraction=0; return; end
    edges=linspace(a,b,max(2,ceil((b-a)*8*f)+1))';
end
power=0;
for k=1:numel(edges)-1
    [tq,wq]=gaussFour(edges(k),edges(k+1));
    power=power+sum(wq.*wave.derivative(tq,m).^2);
end
fraction=power/wave.derivativeEnergy(m);
end

function frac=angularFractions(src,elevationEdges)
% Exact low-degree sphere quadrature for force/moment polynomial patterns.
muEdges=sin(elevationEdges); az=2*pi*(0:7)/8;
mass=zeros(2,numel(muEdges)-1);
for sector=1:size(mass,2)
    [mu,w]=gaussFour(muEdges(sector),muEdges(sector+1));
    [AZ,MU]=meshgrid(az,mu);
    radial=sqrt(max(0,1-MU.^2));
    n=[radial(:).*cos(AZ(:)) radial(:).*sin(AZ(:)) MU(:)];
    weights=repmat(w(:),1,numel(az))*(2*pi/numel(az));
    weights=weights(:);
    if strcmp(src.type,'point_force')
        vector=repmat(src.F(:)',size(n,1),1);
    else
        vector=n*src.M;
    end
    longitudinal=sum(vector.*n,2);
    mass(1,sector)=sum(weights.*longitudinal.^2);
    mass(2,sector)=sum(weights.*max(0,sum(vector.^2,2)-longitudinal.^2));
end

frac=zeros(size(mass));
for beta=1:2
    total=sum(mass(beta,:));
    if total>0, frac(beta,:)=mass(beta,:)/total; end
end
end

function [nodes,weights]=gaussFour(a,b)
% Four-point Gauss-Legendre rule, exact through degree seven.
u=sqrt((3+2*sqrt(6/5))/7);
v=sqrt((3-2*sqrt(6/5))/7);
nodes=(a+b)/2+(b-a)/2*[-u;-v;v;u];
weights=(b-a)/2*[(18-sqrt(30))/36;(18+sqrt(30))/36; ...
    (18+sqrt(30))/36;(18-sqrt(30))/36];
end

function mat=zeroMaterial(vp,vs,rho,f)
mat=MaterialClass; mat.acoustics=false; mat.d=3;
mat.vp=vp; mat.vs=vs; mat.rho=rho; mat.Frequency=f;
mat.CorrelationLength=.01; mat.coefficients_of_variation=[0 0 0];
mat.correlation_coefficients=[0 0 0]; mat.Phi=@(q)ones(size(q));
mat.timeSteps=1; mat.Q=[Inf Inf];
end

function restore(oldPath,oldRng)
path(oldPath); rng(oldRng);
end
