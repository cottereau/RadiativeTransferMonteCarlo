function comparison = RickerPointForceComparison(opts)
%RICKERPOINTFORCECOMPARISON RTE versus exact P/S ballistic point-force cells.
if nargin<1, opts=struct; end
comparison=compareRickerBallisticSource('F',opts);
end
