function src = sourcePointForce(F, vp, vs, rho, opts)
%SOURCEPOINTFORCE  Wave-equation equivalent directional point force source.
%
%   src = sourcePointForce(F, vp, vs, rho, opts) builds a source
%   representation for a single point force vector F = [Fx, Fy, Fz],
%   accounting for elastodynamic radiation patterns and energy partition.
%
%   PHYSICS (Aki & Richards, 2002, Quantitative Seismology, Section 4.1):
%     P-wave far-field radiation amplitude:  A_P(gamma) = (gamma . F)
%     P-wave energy radiation pattern:       R_P(gamma) = (gamma . F)^2
%     S-wave radiation vector:               A_S(gamma) = F - (gamma . F)*gamma
%     S-wave energy radiation pattern:       R_S(gamma) = ||A_S(gamma)||^2
%     S-wave polarization:                   perp = A_S / ||A_S||
%
%     Elastodynamic energy partition (velocity scaling 2*(vP/vS)^3):
%       int R_P dOmega = (4*pi / 3) * ||F||^2
%       int R_S dOmega = (8*pi / 3) * ||F||^2
%       E_S / E_P = 2 * (vP / vS)^3
%
%   INPUTS:
%     F    - 3D force vector [Fx Fy Fz]
%     vp   - P-wave velocity in source medium
%     vs   - S-wave velocity in source medium
%     rho  - Density in source medium (optional, default 1)
%     opts - Struct with optional fields (.position, .lambda, .numberParticles)
%
%   OUTPUT:
%     src - Struct ready for initializeParticle / radiativeTransfer.
%
%   See also SOURCEMOMENTTENSOR, INITIALIZEPARTICLE, RADIATIVETRANSFER.

if nargin < 4 || isempty(rho), rho = 1.0; end
if nargin < 5, opts = struct; end

F_vec = F(:)';
F_norm = norm(F_vec);
if F_norm == 0
    error('Point force vector F cannot be zero.');
end
f_hat = F_vec / F_norm;

if ~isfield(opts,'position') || isempty(opts.position), opts.position = [0 0 0]; end
if ~isfield(opts,'lambda')   || isempty(opts.lambda),   opts.lambda = 1e-3; end
if ~isfield(opts,'numberParticles') || isempty(opts.numberParticles), opts.numberParticles = 1e5; end

% Energy partition
ratio_SP = 2 * (vp / vs)^3;
fP = 1 / (1 + ratio_SP);
fS = ratio_SP / (1 + ratio_SP);

patternP = @(gamma) (sum(gamma .* repmat(f_hat, size(gamma,1), 1), 2)).^2;
patternS = @(gamma) sum((repmat(f_hat, size(gamma,1), 1) - ...
    (sum(gamma .* repmat(f_hat, size(gamma,1), 1), 2)).*gamma).^2, 2);

src = struct('type', 'point_force', ...
             'F', F_vec, ...
             'f_hat', f_hat, ...
             'position', opts.position, ...
             'lambda', opts.lambda, ...
             'numberParticles', opts.numberParticles, ...
             'vp', vp, 'vs', vs, 'rho', rho, ...
             'fP', fP, 'fS', fS, ...
             'ratio_SP', ratio_SP, ...
             'patternP', patternP, ...
             'patternS', patternS);

src.sampleParticles = @(N) sampleForceParticles(N, f_hat, vp, vs, fP, opts);

end

% -------------------------------------------------------------------------
function P = sampleForceParticles(N, f_hat, vp, vs, fP, opts)
NP = round(fP * N);
NS = N - NP;

% P-wave sampling: cos(theta) ~ cos^2(theta) relative to f_hat
% We can sample cos(theta) using inverse CDF or rejection
P_dirs = zeros(NP, 3);
nP_done = 0;
while nP_done < NP
    batch = min(10000, 3*(NP - nP_done));
    c0 = 2*rand(batch,1)-1; ph = 2*pi*rand(batch,1); s0 = sqrt(1-c0.^2);
    cand = [s0.*cos(ph) s0.*sin(ph) c0];
    val = (cand * f_hat').^2;
    acc = rand(batch, 1) < val;
    nAcc = nnz(acc);
    if nAcc > 0
        nTake = min(nAcc, NP - nP_done);
        candAcc = cand(acc, :);
        P_dirs(nP_done+1:nP_done+nTake, :) = candAcc(1:nTake, :);
        nP_done = nP_done + nTake;
    end
end

P_perps = zeros(NP, 3);
if NP > 0
    [~, refIdx] = min(abs(P_dirs), [], 2);
    refVec = zeros(NP, 3);
    refVec(sub2ind(size(refVec), (1:NP)', refIdx)) = 1;
    P_perps = cross(P_dirs, refVec, 2);
    P_perps = P_perps ./ vecnorm(P_perps, 2, 2);
end

% S-wave sampling: sin^2(theta) relative to f_hat
S_dirs = zeros(NS, 3);
S_perps = zeros(NS, 3);
nS_done = 0;
while nS_done < NS
    batch = min(10000, 2*(NS - nS_done));
    c0 = 2*rand(batch,1)-1; ph = 2*pi*rand(batch,1); s0 = sqrt(1-c0.^2);
    cand = [s0.*cos(ph) s0.*sin(ph) c0];
    dot_val = cand * f_hat';
    A_S = repmat(f_hat, batch, 1) - dot_val .* cand;
    val = sum(A_S.^2, 2);
    acc = rand(batch, 1) < val;
    nAcc = nnz(acc);
    if nAcc > 0
        nTake = min(nAcc, NS - nS_done);
        candAcc = cand(acc, :);
        AS_Acc  = A_S(acc, :);
        S_dirs(nS_done+1:nS_done+nTake, :) = candAcc(1:nTake, :);
        normAS = vecnorm(AS_Acc(1:nTake, :), 2, 2);
        S_perps(nS_done+1:nS_done+nTake, :) = AS_Acc(1:nTake, :) ./ max(normAS, eps);
        nS_done = nS_done + nTake;
    end
end

x_pos = repmat(opts.position(:)', N, 1);
if isfield(opts, 'lambda') && opts.lambda > 0
    x_pos = x_pos + opts.lambda * randn(N, 3);
end

P.d = 3;
P.N = N;
P.x = x_pos;
P.dir = [P_dirs; S_dirs];
P.perp = [P_perps; S_perps];
P.p = [true(NP,1); false(NS,1)];
P.t = zeros(N, 1);
P.alive = true(N, 1);
P.dir0 = P.dir;

end
