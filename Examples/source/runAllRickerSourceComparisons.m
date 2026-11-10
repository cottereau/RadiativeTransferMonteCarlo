function results = runAllRickerSourceComparisons(opts)
%RUNALLRICKERSOURCECOMPARISONS Five article-style source mechanisms.
% Each active P/S channel has at most eight particles and three times.
if nargin<1, opts=struct; end
if ~isstruct(opts) || ~isscalar(opts), error('Fonte:ExampleOptions','Options must be a scalar structure.'); end
if ~isfield(opts,'seed'), opts.seed=1801; end
labels={'F','DC','EX','CLVD','GENERAL'};
results=cell(numel(labels),1);
for j=1:numel(labels)
    local=opts; local.seed=opts.seed+j-1;
    if j==1, result=RickerPointForceComparison(local);
    else, result=RickerMomentTensorComparison(labels{j},local); end
    results{j}=result;
    sError=NaN; if ~isempty(result.S), sError=result.S.relativeCellRMSE; end
    fprintf('%s: P=%.9g J, S=%.9g J, S/P=%.9g; cell RMSE P=%.3g, S=%.3g.\n', ...
        result.caseID,result.spectrum.energy(1),result.spectrum.energy(2), ...
        result.expectedRatioSP,result.P.relativeCellRMSE,sError);
end
end
