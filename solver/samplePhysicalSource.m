function P = samplePhysicalSource(N, kind, mechanism, fP, opts, deviator)
% Unit-weight source particles. Mechanisms are scaled to avoid overflow.
% S moment sampling uses only the normalized deviator: isotropic stress
% contributes no transverse amplitude, even for nearly explosive tensors.
validateattributes(N, {'numeric'}, ...
    {'real','finite','scalar','integer','nonnegative','<=',flintmax});
N = double(N);
switch opts.launchMode
    case 'p', NP = N;
    case 's', NP = 0;
    otherwise, NP = unbiasedRound(fP*N);
end
NS = N-NP;
if strcmp(kind,'moment_tensor')
    D = deviator;
    ds = max(abs(D(:)));
    if ds == 0 && strcmp(opts.launchMode,'s')
        error('Fonte:AbsentMode', 'An isotropic moment tensor emits no S energy.');
    end
    if ds > 0, D = D/ds; end
else
    D = mechanism;
end
[pdir, ~] = drawDirections(NP, false, kind, mechanism);
[sdir, sperp] = drawDirections(NS, true, kind, D);
pperp = zeros(NP,3);
if NP > 0
    [~, j] = min(abs(pdir),[],2);
    ref = zeros(NP,3);
    ref(sub2ind(size(ref),(1:NP)',j)) = 1;
    pperp = cross(pdir,ref,2);
    pperp = pperp./sqrt(sum(pperp.^2,2));
end
x = repmat(opts.position,N,1);
if opts.lambda > 0, x = x+opts.lambda*randn(N,3); end
dirs = [pdir;sdir];
P = struct('d',3,'N',N,'x',x,'dir',dirs,'perp',[pperp;sperp], ...
    'p',[true(NP,1);false(NS,1)],'t',zeros(N,1), ...
    'alive',true(N,1),'dir0',dirs);
end

function [dirs, perps] = drawDirections(N, shear, kind, mechanism)
dirs = zeros(N,3); perps = zeros(N,3);
if N == 0, return; end
if strcmp(kind,'point_force')
    bound = sum(mechanism.^2);
else
    bound = norm(mechanism,2)^2;
end
if ~(isfinite(bound) && bound > 0)
    error('Fonte:AbsentMode', 'Cannot sample a zero or invalid radiation pattern.');
end
done = 0; attempts = 0;
% Physical quadratic patterns have acceptance bounded away from zero after
% scaling. Cap attempts to report invalid/numerically broken patterns.
limit = 200*max(N,10);
while done < N
    batch = min([10000, 6*(N-done), limit-attempts]);
    if batch <= 0
        error('Fonte:SamplingFailed', 'Radiation sampler exhausted its candidate budget.');
    end
    attempts = attempts+batch;
    z = 2*rand(batch,1)-1; phi = 2*pi*rand(batch,1);
    r = sqrt(max(0,1-z.^2));
    n = [r.*cos(phi), r.*sin(phi), z];
    if strcmp(kind,'point_force')
        amplitude = repmat(mechanism,batch,1);
    else
        amplitude = n*mechanism;
    end
    normal = sum(amplitude.*n,2);
    if shear
        amplitude = amplitude-normal.*n;
        value = sum(amplitude.^2,2);
    else
        value = normal.^2;
    end
    selected = find(value > 0 & rand(batch,1) < min(1,value/bound));
    take = min(numel(selected),N-done);
    selected = selected(1:take);
    dirs(done+(1:take),:) = n(selected,:);
    if shear
        perps(done+(1:take),:) = amplitude(selected,:)./sqrt(value(selected));
    end
    done = done+take;
end
end
