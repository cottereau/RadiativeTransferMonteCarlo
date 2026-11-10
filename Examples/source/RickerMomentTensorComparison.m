function comparison = RickerMomentTensorComparison(mechanism,opts)
%RICKERMOMENTTENSORCOMPARISON DC, EX, CLVD or GENERAL Ricker comparison.
if nargin<1, mechanism='DC'; end
if nargin<2, opts=struct; end
mechanism=upper(char(mechanism));
if ~any(strcmp(mechanism,{'DC','EX','CLVD','GENERAL'}))
    error('Fonte:ExampleMechanism','Use DC, EX, CLVD or GENERAL.');
end
comparison=compareRickerBallisticSource(mechanism,opts);
end
