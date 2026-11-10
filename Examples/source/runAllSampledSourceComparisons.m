function results = runAllSampledSourceComparisons(time,values,opts)
%RUNALLSAMPLEDSOURCECOMPARISONS Arbitrary tabulated pulse for all sources.
% If time/values are omitted, use a short illustrative nonanalytic series.
if nargin<1, time=[-.03 -.02 -.01 0 .011 .021 .032]; end
if nargin<2, values=[0 .2 -.6 1 -.2 .1 0]; end
if nargin<3, opts=struct; end
if ~isstruct(opts) || ~isscalar(opts), error('Fonte:ExampleOptions','Options must be a scalar structure.'); end
if ~isfield(opts,'seed'), opts.seed=1901; end
labels={'F','DC','EX','CLVD','GENERAL'};
results=cell(numel(labels),1);
for j=1:numel(labels)
    local=opts; local.seed=opts.seed+j-1;
    if j==1, result=SampledPointForceComparison(time,values,local);
    else, result=SampledMomentTensorComparison(labels{j},time,values,local); end
    results{j}=result;
    sError=NaN; if ~isempty(result.S), sError=result.S.relativeCellRMSE; end
    fprintf('%s: P=%.9g J, S=%.9g J, S/P=%.9g; cell RMSE P=%.3g, S=%.3g.\n', ...
        result.caseID,result.spectrum.energy(1),result.spectrum.energy(2), ...
        result.expectedRatioSP,result.P.relativeCellRMSE,sError);
end
end
