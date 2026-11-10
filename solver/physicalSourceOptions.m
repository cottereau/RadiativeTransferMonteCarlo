function opts = physicalSourceOptions(vp, vs, rho, opts)
% Shared contract for 3D bulk sources.
% vs = [] denotes an acoustic source; vp is then the sound speed.
validateattributes(vp, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(rho, {'numeric'}, {'real','finite','scalar','positive'});
if ~isempty(vs)
    validateattributes(vs, {'numeric'}, {'real','finite','scalar','positive'});
    if double(vp)/double(vs) <= sqrt(4/3)
        error('Fonte:InvalidMedium', 'Positive bulk modulus requires vp/vs > sqrt(4/3).');
    end
end
if ~isstruct(opts) || ~isscalar(opts)
    error('Fonte:InvalidOptions', 'Source options must be a scalar structure.');
end
defaults = struct('position',[0 0 0], 'lambda',1e-3, ...
    'numberParticles',1e5, 'launchMode','mixed');
names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
    if ~isfield(opts,name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
validateattributes(opts.position, {'numeric'}, {'real','finite','vector','numel',3});
validateattributes(opts.lambda, {'numeric'}, {'real','finite','scalar','nonnegative'});
validateattributes(opts.numberParticles, {'numeric'}, ...
    {'real','finite','scalar','integer','nonnegative','<=',flintmax});
mode = opts.launchMode;
if isstring(mode) && isscalar(mode), mode = char(mode); end
if ~ischar(mode) || ~isrow(mode) || ~any(strcmpi(mode, {'P','S','mixed'}))
    error('Fonte:InvalidLaunchMode', 'launchMode must be P, S, or mixed.');
end
opts.launchMode = lower(mode);
opts.position = double(opts.position(:)');
opts.lambda = double(opts.lambda);
opts.numberParticles = double(opts.numberParticles);
if isfield(opts,'halfspace') && ~isempty(opts.halfspace)
    error('Fonte:UnsupportedHalfspace', ...
        'Bulk physical sources do not implement half-space radiation.');
end
end
