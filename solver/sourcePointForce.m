function src = sourcePointForce(F, vp, vs, rho, opts)
%SOURCEPOINTFORCE  Normalized 3D bulk point-force radiation and P/S fractions.
%   F is a real nonzero 3-vector in N. vp, vs, rho use consistent SI units.
%   opts: position, lambda (Cartesian jitter std), numberParticles, launchMode.
%   launchMode = 'mixed' (default), 'P', or 'S'; fP/fS always describe the
%   physical mixed source. Conditional launches preserve the angular law.
%   This helper does not yet assign absolute waveform-dependent energy.
if nargin < 4 || isempty(rho), rho = 1; end
if nargin < 5, opts = struct; end
opts = physicalSourceOptions(vp,vs,rho,opts);
validateattributes(F, {'numeric'}, {'real','finite','vector','numel',3});
F = double(F(:)');
scale = max(abs(F));
if scale == 0, error('Fonte:ZeroMechanism', 'Point force must be nonzero.'); end
fhat = F/scale;
fhat = fhat/norm(fhat);
% Scale by vs^3: avoids overflow from (vp/vs)^3 in the fractions.
p = (double(vs)/double(vp))^3;
fP = p/(p+2); fS = 2/(p+2);
src = struct('type','point_force','F',F,'f_hat',fhat, ...
    'position',opts.position,'lambda',opts.lambda, ...
    'numberParticles',opts.numberParticles,'launchMode',opts.launchMode, ...
    'vp',double(vp),'vs',double(vs),'rho',double(rho), ...
    'fP',fP,'fS',fS,'ratio_SP',fS/fP, ...
    'patternP',@(n) (n*fhat').^2, ...
    'patternS',@(n) sum((fhat-(n*fhat').*n).^2,2));
src.sampleParticles = @(N) samplePhysicalSource(N,'point_force',fhat,fP,opts);
end
