function src = sourceMomentTensor(M, vp, vs, rho, opts)
%SOURCEMOMENTTENSOR  Wave-equation equivalent seismic moment tensor source.
%
%   src = sourceMomentTensor(M, vp, vs, rho, opts) builds a source
%   representation for an arbitrary 3x3 symmetric seismic moment tensor M
%   (e.g. double-couple fault, explosive/implosive, CLVD, or tensile crack),
%   accounting for elastodynamic radiation patterns and energy partition.
%
%   PHYSICS (Aki & Richards, 2002, Quantitative Seismology, Ch. 4):
%     P-wave far-field radiation amplitude:  A_P(gamma) = gamma' * M * gamma
%     P-wave energy radiation pattern:       R_P(gamma) = (gamma' * M * gamma)^2
%     S-wave radiation vector:               A_S(gamma) = M*gamma - (gamma'*M*gamma)*gamma
%     S-wave energy radiation pattern:       R_S(gamma) = ||A_S(gamma)||^2
%     S-wave polarization:                   perp = A_S / ||A_S||
%
%     Elastodynamic energy partition (velocity scaling (vP/vS)^5):
%       E_P proportional to (1 / (rho * vP^5)) * int R_P dOmega
%       E_S proportional to (1 / (rho * vS^5)) * int R_S dOmega
%
%   INPUTS:
%     M    - 3x3 symmetric matrix (or [strike, dip, rake, M0] in degrees/N.m)
%     vp   - P-wave velocity in source medium
%     vs   - S-wave velocity in source medium
%     rho  - Density in source medium (optional, default 1)
%     opts - Struct with optional fields:
%              .position  [x y z] (default [0 0 0])
%              .lambda    Gaussian source radius (default 1e-3)
%              .numberParticles (default 1e5)
%
%   OUTPUT:
%     src - Struct ready for initializeParticle / radiativeTransfer, with:
%             .type = 'point'
%             .momentTensor = M
%             .fP, .fS (energy fractions)
%             .patternP, .patternS (radiation pattern functions)
%             .sampleInitial (function handle generating non-isotropic P)
%
%   See also SOURCEPOINTFORCE, INITIALIZEPARTICLE, RADIATIVETRANSFER.

if nargin < 4 || isempty(rho), rho = 1.0; end
if nargin < 5, opts = struct; end

% Handle fault orientation input [strike, dip, rake, M0]
if numel(M) == 4
    phi = M(1)*pi/180; delta = M(2)*pi/180; lambda = M(3)*pi/180; M0 = M(4);
    % Fault normal and slip vector (Aki & Richards convention)
    n = [-sin(delta)*sin(phi), -sin(delta)*cos(phi), -cos(delta)];
    s = [cos(lambda)*cos(phi) + cos(delta)*sin(lambda)*sin(phi), ...
        -cos(lambda)*sin(phi) + cos(delta)*sin(lambda)*cos(phi), ...
        -sin(delta)*sin(lambda)];
    M = M0 * (n(:)*s(:).' + s(:)*n(:).');
end

% Ensure symmetry
M = (M + M.') / 2;

% Default options
if ~isfield(opts,'position') || isempty(opts.position), opts.position = [0 0 0]; end
if ~isfield(opts,'lambda')   || isempty(opts.lambda),   opts.lambda = 1e-3; end
if ~isfield(opts,'numberParticles') || isempty(opts.numberParticles), opts.numberParticles = 1e5; end

% Radiation pattern function handles
patternP = @(gamma) (sum((gamma * M) .* gamma, 2)).^2;
patternS = @(gamma) sum((gamma * M - (sum((gamma * M).*gamma, 2)).*gamma).^2, 2);

% Integrate total P and S energy over unit sphere by numerical quadrature
Nquad = 50000;
c0 = 2*rand(Nquad,1)-1; ph = 2*pi*rand(Nquad,1); s0 = sqrt(1-c0.^2);
dirs = [s0.*cos(ph) s0.*sin(ph) c0];

int_RP = mean(patternP(dirs)); % int R_P dOmega / 4*pi
int_RS = mean(patternS(dirs)); % int R_S dOmega / 4*pi

rad_energy_P = int_RP / (rho * vp^5);
rad_energy_S = int_RS / (rho * vs^5);
if rad_energy_S < 1e-12 * rad_energy_P
    rad_energy_S = 0;
end
if rad_energy_P < 1e-12 * rad_energy_S
    rad_energy_P = 0;
end
tot_energy = rad_energy_P + rad_energy_S;

if tot_energy > 0
    fP = rad_energy_P / tot_energy;
    fS = rad_energy_S / tot_energy;
else
    fP = 1; fS = 0;
end

src = struct('type', 'moment_tensor', ...
             'M', M, ...
             'position', opts.position, ...
             'lambda', opts.lambda, ...
             'numberParticles', opts.numberParticles, ...
             'vp', vp, 'vs', vs, 'rho', rho, ...
             'fP', fP, 'fS', fS, ...
             'patternP', patternP, ...
             'patternS', patternS);

% Generator function for sampling N particles according to M
src.sampleParticles = @(N) sampleMomentParticles(N, M, vp, vs, fP, opts);

end

% -------------------------------------------------------------------------
function P = sampleMomentParticles(N, M, vp, vs, fP, opts)
% Generate N particles distributed between P and S waves according to the
% moment tensor radiation patterns and energy partition fP, fS.

NP = round(fP * N);
NS = N - NP;

% Sample P-wave particles using rejection sampling on patternP
P_dirs = zeros(NP, 3);
P_perps = zeros(NP, 3);
maxP = norm(M, 2)^2; % upper bound on patternP

nP_done = 0;
while nP_done < NP
    batch = min(10000, 3*(NP - nP_done));
    c0 = 2*rand(batch,1)-1; ph = 2*pi*rand(batch,1); s0 = sqrt(1-c0.^2);
    cand = [s0.*cos(ph) s0.*sin(ph) c0];
    val = (sum((cand * M) .* cand, 2)).^2;
    acc = rand(batch, 1) < (val / maxP);
    nAcc = nnz(acc);
    if nAcc > 0
        nTake = min(nAcc, NP - nP_done);
        candAcc = cand(acc, :);
        P_dirs(nP_done+1:nP_done+nTake, :) = candAcc(1:nTake, :);
        nP_done = nP_done + nTake;
    end
end

% Orthonormal perp vectors for P (arbitrary in tangent plane)
if NP > 0
    [~, refIdx] = min(abs(P_dirs), [], 2);
    refVec = zeros(NP, 3);
    refVec(sub2ind(size(refVec), (1:NP)', refIdx)) = 1;
    P_perps = cross(P_dirs, refVec, 2);
    P_perps = P_perps ./ vecnorm(P_perps, 2, 2);
end

% Sample S-wave particles using rejection sampling on patternS
S_dirs = zeros(NS, 3);
S_perps = zeros(NS, 3);
maxS = norm(M, 2)^2;

nS_done = 0;
while nS_done < NS
    batch = min(10000, 3*(NS - nS_done));
    c0 = 2*rand(batch,1)-1; ph = 2*pi*rand(batch,1); s0 = sqrt(1-c0.^2);
    cand = [s0.*cos(ph) s0.*sin(ph) c0];
    A_S = cand * M - (sum((cand * M) .* cand, 2)) .* cand;
    val = sum(A_S.^2, 2);
    acc = rand(batch, 1) < (val / maxS);
    nAcc = nnz(acc);
    if nAcc > 0
        nTake = min(nAcc, NS - nS_done);
        candAcc = cand(acc, :);
        AS_Acc  = A_S(acc, :);
        S_dirs(nS_done+1:nS_done+nTake, :) = candAcc(1:nTake, :);
        % Polarize S particle along A_S
        normAS = vecnorm(AS_Acc(1:nTake, :), 2, 2);
        S_perps(nS_done+1:nS_done+nTake, :) = AS_Acc(1:nTake, :) ./ max(normAS, eps);
        nS_done = nS_done + nTake;
    end
end

% Combine P and S populations
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
