function counts = sourceParticleCounts(sources,N)
% Globally allocate unit-weight particles; no packet-dependent rounding.
if ~iscell(sources) || isempty(sources)
    error('Fonte:InvalidSource','A nonempty cell array of sources is required.');
end
f = zeros(numel(sources),1);
for j = 1:numel(sources)
    s = sources{j};
    if ~isstruct(s) || ~isscalar(s) || ~isfield(s,'fraction')
        error('Fonte:MissingFraction','Every source must have a scalar fraction.');
    end
    validateattributes(s.fraction,{'numeric'},{'real','finite','scalar','nonnegative'});
    f(j) = s.fraction;
end
if ~any(f), error('Fonte:ZeroFractions','At least one fraction must be positive.'); end
if abs(sum(f)-1)>1e-6
    warning('Fonte:NormalizedFractions','Source fractions do not sum to one; normalizing.');
end
f = f/max(f); f = f/sum(f);
expected = N*f; counts = floor(expected);
[~,order] = sortrows([-(expected-counts),(1:numel(f))'],[1 2]);
remaining = N-sum(counts);
counts(order(1:remaining)) = counts(order(1:remaining))+1;
end
