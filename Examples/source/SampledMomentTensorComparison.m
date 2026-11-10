function comparison = SampledMomentTensorComparison(mechanism,time,values,opts)
%SAMPLEDMOMENTTENSORCOMPARISON DC, EX, CLVD or GENERAL sampled pulse.
% time in seconds, values dimensionless; moment amplitude is opts.amplitude Nm.
if nargin<4, opts=struct; end
mechanism=upper(char(mechanism));
if ~any(strcmp(mechanism,{'DC','EX','CLVD','GENERAL'}))
    error('Fonte:ExampleMechanism','Use DC, EX, CLVD or GENERAL.');
end
opts.sampledTime=time; opts.sampledValues=values;
comparison=compareRickerBallisticSource(mechanism,opts);
end
