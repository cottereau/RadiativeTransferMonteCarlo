function src = sourceMomentTensor(M, vp, vs, rho, opts)
%SOURCEMOMENTTENSOR  Exact normalized radiation for a real symmetric 3D moment.
%   M: 3x3 tensor in N m, or [strike dip rake M0] (angles in degrees) using
%   Aki & Richards (2002) Box 4.4: x = north, y = east, z = down.
%   opts: position, lambda (Cartesian jitter std), numberParticles, launchMode.
%   launchMode = 'mixed' (default), 'P', or 'S'. An isotropic moment has no S
%   conditional distribution. fP/fS retain physical mixed-source fractions.
%   No random quadrature or amplitude threshold is used in construction.
%   patternP/S include M's amplitude; sampling uses scaled mechanisms.
%
%   References:
%     - Aki, K. & Richards, P.G. (2002), Quantitative Seismology, 2nd ed., University
%       Science Books (no DOI; ISBN 978-0-935702-96-5)
if nargin < 4 || isempty(rho), rho = 1; end
if nargin < 5, opts = struct; end
opts = physicalSourceOptions(vp,vs,rho,opts);
validateattributes(M, {'numeric'}, {'real','finite','nonempty'});
M = double(M);
if isvector(M) && numel(M) == 4
    phi = M(1)*pi/180; delta = M(2)*pi/180;
    rake = M(3)*pi/180; M0 = M(4);
    % Fault normal nu and slip d (Aki & Richards 2002, Box 4.4).
    nu = [-sin(delta)*sin(phi), sin(delta)*cos(phi), -cos(delta)];
    d = [cos(rake)*cos(phi)+cos(delta)*sin(rake)*sin(phi), ...
        cos(rake)*sin(phi)-cos(delta)*sin(rake)*cos(phi), ...
        -sin(rake)*sin(delta)];
    M = M0*(nu(:)*d(:).'+d(:)*nu(:).');
end
if ~isequal(size(M),[3 3]) || any(~isfinite(M(:)))
    error('Fonte:InvalidMoment', 'Moment must be a finite 3x3 matrix or four fault parameters.');
end
scale = max(abs(M(:)));
if scale == 0, error('Fonte:ZeroMechanism', 'Moment tensor must be nonzero.'); end
Q = M/scale;
if norm(Q-Q.','fro') > 1e-12*norm(Q,'fro')
    error('Fonte:NonSymmetricMoment', 'Moment tensor must be symmetric.');
end
M = M/2+M.'/2;
Q = Q/2+Q.'/2;
% Pairwise differences retain tiny diagonal anisotropy that subtracting
% a rounded trace/3 can lose. Same-sign subtraction cannot overflow.
D = sourceMomentDeviator(M,scale);
CP = trace(Q)^2+2*sum(Q(:).^2);
CS = 3*sum(D(:).^2);
p = CP*(double(vs)/double(vp))^5;
fP = p/(p+CS); fS = CS/(p+CS);
if ~any(D(:))
    fP = 1; fS = 0;
    if strcmp(opts.launchMode,'s')
        error('Fonte:AbsentMode', 'An isotropic moment tensor emits no S energy.');
    end
end
% Removing the isotropic part before transverse projection avoids cancellation.
Dphysical = D*scale;
% Use normalized D for sampling, and the physical tensor for public patterns.
src = struct('type','moment_tensor','M',M, ...
    'position',opts.position,'lambda',opts.lambda, ...
    'numberParticles',opts.numberParticles,'launchMode',opts.launchMode, ...
    'vp',double(vp),'vs',double(vs),'rho',double(rho), ...
    'fP',fP,'fS',fS,'ratio_SP',fS/fP, ...
    'patternP',@(n) sum((n*M).*n,2).^2, ...
    'patternS',@(n) sum((n*Dphysical-sum((n*Dphysical).*n,2).*n).^2,2));
src.sampleParticles = @(N) samplePhysicalSource(N,'moment_tensor',Q,fP,opts,D);
end
