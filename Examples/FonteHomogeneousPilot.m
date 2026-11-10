function pilot = FonteHomogeneousPilot
%FONTEHOMOGENEOUSPILOT Tiny full-spectrum ballistic Fonte source check.
% Manuscript's proposed homogeneous pre-SEM benchmark: F-R, DC-R, EX-R.
% The numerical background below is illustrative, not a manuscript fit.
% In zero scattering, conditional transport is frequency-independent, so
% the exact integrated source energies weight one P/S response per case.
root=fileparts(fileparts(mfilename('fullpath')));
if isempty(which('sourceEnergySpectrum')), addpath(root,'-begin'); AddRTMCLib(root,false); end
geometry=struct('dimension',3,'frame','spherical');
observation=struct('x',[0 25 50 100],'y',[-pi pi], ...
    'z',[-pi/2 pi/2],'directions',[0 pi],'time',[0 .01 .02]);
solverOptions=struct('eventTime',0,'packetSize',4, ...
    'scatteringAngleNodes',16,'polarized',true, ...
    'maxPolarizedParticles',4);
ids={'F-R','DC-R','EX-R'};
pilot=struct('caseID',{},'spectrumEnergy',{},'expectedRatioSP',{}, ...
    'energy',{},'time',{},'source',{},'waveform',{},'result',{});
for j=1:numel(ids)
    if j==1, example=FontePointForceBroadband('R');
    elseif j==2, example=FonteMomentTensorBroadband('DC','R');
    else, example=FonteMomentTensorBroadband('EX','R'); end
    spec=example.spectrum;
    responses=cell(1,2);
    for beta=1:2
        if spec.energy(beta)==0, continue; end
        src=spec.source; src.launchMode='P';
        if beta==2, src.launchMode='S'; end
        responses{beta}=radiativeTransfer(geometry,src, ...
            homogeneousMaterial(10),observation,solverOptions);
    end
    physical=struct('energyP',spec.energy(1),'energyS',spec.energy(2), ...
        'eventTime',spec.eventTime);
    result=composeSourceBands(physical,responses);
    total=sum(spec.energy);
    if any(abs(sum(result.energy,2)-total)>1e-11*total)
        error('Fonte:PilotEnergy','All-space energy did not match the physical source.');
    end
    ratio=spec.energy(2)/spec.energy(1);
    pilot(j)=struct('caseID',ids{j},'spectrumEnergy',spec.energy, ...
        'expectedRatioSP',ratio,'energy',result.energy,'time',result.t, ...
        'source',spec.source,'waveform',example.waveform,'result',result);
    fprintf('%s: total energy %.9g J; S/P %.9g; %d active channels.\n', ...
        ids{j},total,ratio,nnz(spec.energy));
end
end

function mat=homogeneousMaterial(f)
mat=MaterialClass; mat.acoustics=false; mat.d=3;
mat.vp=4000; mat.vs=2000; mat.rho=2500; mat.Frequency=f;
mat.CorrelationLength=.01; mat.coefficients_of_variation=[0 0 0];
mat.correlation_coefficients=[0 0 0]; mat.Phi=@(q)ones(size(q));
mat.timeSteps=1; mat.Q=[Inf Inf];
end
