function example = fonteBenchmarkCase(mechanism,waveKind,opts)
%FONTEBENCHMARKCASE Reproducible source configuration for Fonte cases.
% Configuration is read-only by default. Set opts.run=true together with
% a materialFactory and explicitly budgeted bands for transport.
if nargin<3, opts=struct; end
if ~isstruct(opts) || ~isscalar(opts), error('Fonte:ExampleOptions','Options must be a scalar structure.'); end
defaults=struct('vp',4000,'vs',2000,'rho',2500,'amplitude',1, ...
    'numberParticles',4,'tensor',[],'run',false,'materialFactory',[], ...
    'bands',[],'geometry',[],'observation',[],'solverOptions',[], ...
    'maxRuns',8,'maxTotalParticles',100,'seed',1706);
names=fieldnames(defaults);
if any(~ismember(fieldnames(opts),names)), error('Fonte:ExampleOptions','Unknown example option.'); end
for j=1:numel(names)
    if ~isfield(opts,names{j}), opts.(names{j})=defaults.(names{j}); end
end
validateattributes(opts.amplitude,{'numeric'},{'real','finite','scalar','nonzero'});
validateattributes(opts.numberParticles,{'numeric'},{'real','finite','scalar','integer','positive'});
mechanism=upper(char(mechanism)); waveKind=lower(char(waveKind));
switch waveKind
    case {'r','ricker'}
        wave=sourceTimeFunction('ricker',struct('f0',10,'delay',0)); code='R';
    case {'o','ormsby'}
        wave=sourceTimeFunction('ormsby',struct('corners',[4 8 12 16],'delay',0)); code='O';
    case {'t','toneburst'}
        wave=sourceTimeFunction('toneburst',struct('fc',10,'cycles',3,'delay',0)); code='T';
    otherwise
        error('Fonte:ExampleWaveform','Use R, O or T.');
end
srcopts=struct('numberParticles',opts.numberParticles,'lambda',0);
switch mechanism
    case 'F'
        source=sourcePointForce(opts.amplitude*[0 0 1],opts.vp,opts.vs,opts.rho,srcopts);
    case 'DC'
        M=[0 1 0;1 0 0;0 0 0];
        source=sourceMomentTensor(opts.amplitude*M,opts.vp,opts.vs,opts.rho,srcopts);
    case 'EX'
        source=sourceMomentTensor(opts.amplitude*eye(3),opts.vp,opts.vs,opts.rho,srcopts);
    case 'CLVD'
        source=sourceMomentTensor(opts.amplitude*diag([1 1 -2]),opts.vp,opts.vs,opts.rho,srcopts);
    case 'GENERAL'
        if isempty(opts.tensor), opts.tensor=[1 .2 0;.2 -.3 .1;0 .1 .5]; end
        source=sourceMomentTensor(opts.amplitude*opts.tensor,opts.vp,opts.vs,opts.rho,srcopts);
    otherwise
        error('Fonte:ExampleMechanism','Use F, DC, EX, CLVD or GENERAL.');
end
spectrum=sourceEnergySpectrum(source,wave);
example=struct('caseID',[mechanism '-' code],'source',source, ...
    'waveform',wave,'spectrum',spectrum,'units','SI', ...
    'description','3D homogeneous bulk point source; source spectrum in physical joules');
if ~opts.run, return; end
if ~isa(opts.materialFactory,'function_handle') || isempty(opts.bands)
    error('Fonte:ExampleRun','Running requires a fresh materialFactory and explicit bands.');
end
if isempty(opts.geometry), opts.geometry=struct('dimension',3,'frame','spherical'); end
if isempty(opts.observation)
    opts.observation=struct('x',[0 100],'y',[-pi pi], ...
        'z',[-pi/2 pi/2],'directions',[0 pi], ...
        'time',wave.eventTime(spectrum.derivativeOrder)+[0 .01 .02]);
end
if isempty(opts.solverOptions)
    opts.solverOptions=struct('packetSize',opts.numberParticles, ...
        'scatteringAngleNodes',16,'polarized',true);
end
runopts=struct('bands',opts.bands,'solverOptions',opts.solverOptions, ...
    'maxRuns',opts.maxRuns,'maxTotalParticles',opts.maxTotalParticles, ...
    'seed',opts.seed);
example.result=radiativeTransferBroadband(opts.geometry,spectrum, ...
    opts.materialFactory,opts.observation,runopts);
example.relaxation=fonteModalRelaxation(example.result.responseBands, ...
    example.result.materialMetadata,example.result.t);
end
