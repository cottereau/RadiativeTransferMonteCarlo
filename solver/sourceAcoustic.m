function src = sourceAcoustic(Q, F, c, rho, opts)
%SOURCEACOUSTIC  3D acoustic point source: volume injection and/or force.
%   q = Q s(t) delta(x-x0) (Q in m^3/s) and f = F s(t) delta(x-x0) (F in N)
%   in beta p_t + div v = q, rho v_t + grad p = f. Radiation per unit solid
%   angle is proportional to (rho*c*Q + n.F)^2: monopole, dipole, or their
%   cardioid when both act with the same s(t). One mode only (fP = 1).
%   c is the sound speed (material.v). opts: position, lambda (Cartesian
%   jitter std), numberParticles, launchMode ('mixed' or 'P').
if nargin < 5, opts = struct; end
opts = physicalSourceOptions(c,[],rho,opts);
if strcmp(opts.launchMode,'s')
    error('Fonte:AbsentMode', 'An acoustic source emits no S energy.');
end
validateattributes(Q, {'numeric'}, {'real','finite','scalar'});
if isempty(F), F = [0 0 0]; end
validateattributes(F, {'numeric'}, {'real','finite','vector','numel',3});
Q = double(Q); F = double(F(:)'); c = double(c); rho = double(rho);
a = rho*c*Q;                       % monopole amplitude in force units (N)
scale = max(abs([a F]));
if scale == 0, error('Fonte:ZeroMechanism', 'Acoustic source must be nonzero.'); end
src = struct('type','acoustic','Q',Q,'F',F,'c',c,'rho',rho, ...
    'position',opts.position,'lambda',opts.lambda, ...
    'numberParticles',opts.numberParticles,'launchMode',opts.launchMode, ...
    'fP',1,'fS',0, ...
    'pattern',@(n) (a+n*F').^2/rho);
src.sampleParticles = @(N) sampleAcoustic(N,a/scale,F/scale,opts);
end

function P = sampleAcoustic(N, a, f, opts)
% Rejection from the uniform sphere with bound (|a|+|f|)^2 >= (a+n.f)^2.
validateattributes(N, {'numeric'}, ...
    {'real','finite','scalar','integer','nonnegative','<=',flintmax});
N = double(N); dirs = zeros(N,3); done = 0;
bound = (abs(a)+norm(f))^2;
while done < N
    batch = min(10000, 6*(N-done));
    z = 2*rand(batch,1)-1; phi = 2*pi*rand(batch,1); r = sqrt(max(0,1-z.^2));
    n = [r.*cos(phi), r.*sin(phi), z];
    keep = find(rand(batch,1)*bound < (a+n*f').^2);
    take = min(numel(keep),N-done);
    dirs(done+(1:take),:) = n(keep(1:take),:);
    done = done+take;
end
% Any unit vector orthogonal to dir; acoustic transport does not use it.
[~,j] = min(abs(dirs),[],2);
ref = zeros(N,3); ref(sub2ind(size(ref),(1:N)',j)) = 1;
perp = cross(dirs,ref,2); perp = perp./sqrt(sum(perp.^2,2));
x = repmat(opts.position,N,1);
if opts.lambda > 0, x = x+opts.lambda*randn(N,3); end
P = struct('d',3,'N',N,'x',x,'dir',dirs,'perp',perp,'p',true(N,1), ...
    't',zeros(N,1),'alive',true(N,1),'dir0',dirs);
end
