function study = FonteMultifrequencyVarianceCase
%FONTEMULTIFREQUENCYVARIANCECASE Tiny two-band RTE variance case study.
% F-R source; upper-hemisphere cell has exact P/S particle probability 1/2.
% Two independent batches, four particles/mode/band. No random scattering.
% Contrasts extra particles with a fixed total particle budget.
oldPath=path; oldRng=rng;
cleanup=onCleanup(@()restore(oldPath,oldRng)); %#ok<NASGU>
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root,'-begin'); AddRTMCLib(root,false);
N=4; B=2; probability=.5;
source=sourcePointForce([0 0 1],4000,2000,2500, ...
    struct('numberParticles',N,'lambda',0));
wave=sourceTimeFunction('ricker',struct('f0',10,'delay',0));
spectrum=sourceEnergySpectrum(source,wave);
omega=2*pi*[8;16]; weightOmega=2*pi*[4;4];
p=weightOmega.*spectrum.P(omega); s=weightOmega.*spectrum.S(omega);
bands=struct('omega',omega,'frequencyHz',omega/(2*pi), ...
    'weightOmega',weightOmega,'energyP',p,'energyS',s, ...
    'eventTime',spectrum.eventTime,'retainedTotalEnergy',sum(p)+sum(s));
energy=[p s]; totalByMode=sum(energy,1);
geometry=struct('dimension',3,'frame','spherical');
observation=struct('x',[0 100],'y',[-pi pi], ...
    'z',[-pi/2 0 pi/2],'directions',[0 pi],'time',.01);
solverOptions=struct('packetSize',N,'scatteringAngleNodes',16, ...
    'polarized',true, ...
    'maxPolarizedParticles',N);

% Frequency-independent ballistic response: repeat that same unit problem
% twice, then multiply the average by the total two-node physical energy.
unitCell=zeros(B,2);
for batch=1:B
    for beta=1:2
        conditional=source; conditional.launchMode='P';
        if beta==2, conditional.launchMode='S'; end
        rng(3000+(batch-1)*2+beta,'twister');
        r=radiativeTransfer(geometry,conditional,zeroMaterial(12), ...
            observation,solverOptions);
        unitCell(batch,beta)=r.energyDensity(1,2,1,beta)*r.cellMeasure(1,2);
    end
end
singleMean=mean(unitCell,1).*totalByMode;
singleObservedSE=std(unitCell,0,1)/sqrt(B).*totalByMode;
singleCov=diag(repmat(probability*(1-probability)/(N*B),2,1));
singleBudget=fonteCompositionVariance(totalByMode,singleCov);

driver=struct('bands',bands,'solverOptions',solverOptions, ...
    'batches',B,'seed',3101,'maxRuns',8,'maxTotalParticles',32);
two=radiativeTransferBroadband(geometry,spectrum,@zeroMaterial, ...
    observation,driver);
twoMean=[two.energyDensity(1,2,1,1) two.energyDensity(1,2,1,2)]*two.cellMeasure(1,2);
twoObservedSE=[two.mcStandardError.energyDensity(1,2,1,1) ...
    two.mcStandardError.energyDensity(1,2,1,2)]*two.cellMeasure(1,2);
perRunVariance=probability*(1-probability)/(N*B);
twoBudget=fonteCompositionVariance(energy,diag(repmat(perRunVariance,4,1)));

% Same TOTAL particles as the repeated single-response baseline: two per
% frequency and mode per batch, rather than four at a single frequency.
sourceHalf=source; sourceHalf.numberParticles=N/2;
spectrumHalf=sourceEnergySpectrum(sourceHalf,wave);
fixedSolver=solverOptions; fixedSolver.packetSize=N/2;
fixedSolver.maxPolarizedParticles=N/2;
fixedDriver=driver; fixedDriver.solverOptions=fixedSolver;
fixedDriver.seed=3201; fixedDriver.maxTotalParticles=16;
fixed=radiativeTransferBroadband(geometry,spectrumHalf,@zeroMaterial, ...
    observation,fixedDriver);
fixedMean=[fixed.energyDensity(1,2,1,1) fixed.energyDensity(1,2,1,2)]*fixed.cellMeasure(1,2);
fixedObservedSE=[fixed.mcStandardError.energyDensity(1,2,1,1) ...
    fixed.mcStandardError.energyDensity(1,2,1,2)]*fixed.cellMeasure(1,2);
fixedBudget=fonteCompositionVariance(energy, ...
    diag(repmat(probability*(1-probability)/((N/2)*B),4,1)));

study=struct('bands',bands,'source',source,'waveform',wave, ...
    'upperHemisphereProbability',probability, ...
    'expectedUpperEnergyJ',probability*totalByMode, ...
    'singleRepeated',struct('meanJ',singleMean,'observedSEJ',singleObservedSE, ...
        'exactVariance',singleBudget,'particles',2*B*N), ...
    'twoBandsExtraWork',struct('meanJ',twoMean,'observedSEJ',twoObservedSE, ...
        'exactVariance',twoBudget,'particles',two.particleCount), ...
    'twoBandsFixedWork',struct('meanJ',fixedMean,'observedSEJ',fixedObservedSE, ...
        'exactVariance',fixedBudget,'particles',fixed.particleCount), ...
    'extraWorkVarianceRatio',twoBudget.variance/singleBudget.variance, ...
    'fixedWorkVarianceRatio',fixedBudget.variance/singleBudget.variance, ...
    'effectiveParticlesPerModeExtra',totalByMode.^2./ ...
        (sum(energy.^2,1)/(N*B)), ...
    'effectiveParticlesPerModeFixed',totalByMode.^2./ ...
        (sum(energy.^2,1)/((N/2)*B)), ...
    'note','Exact variances apply to this ballistic hemisphere indicator; two-batch empirical SE is illustrative only.');
fprintf('Variance case: two bands / repeated same case = %.6f with twice the particles.\n', ...
    study.extraWorkVarianceRatio);
fprintf('Variance case: two bands / repeated same case = %.6f at fixed total particles.\n', ...
    study.fixedWorkVarianceRatio);
end

function mat=zeroMaterial(f)
mat=MaterialClass; mat.acoustics=false; mat.d=3;
mat.vp=4000; mat.vs=2000; mat.rho=2500; mat.Frequency=f;
mat.CorrelationLength=.01; mat.coefficients_of_variation=[0 0 0];
mat.correlation_coefficients=[0 0 0]; mat.Phi=@(q)ones(size(q));
mat.timeSteps=1; mat.Q=[Inf Inf];
end

function restore(oldPath,oldRng)
path(oldPath); rng(oldRng);
end
