function comparison = SampledPointForceComparison(time,values,opts)
%SAMPLEDPOINTFORCECOMPARISON Arbitrary measured/tabulated force pulse.
% time in seconds, values dimensionless; force amplitude is opts.amplitude N.
if nargin<3, opts=struct; end
opts.sampledTime=time; opts.sampledValues=values;
comparison=compareRickerBallisticSource('F',opts);
end
