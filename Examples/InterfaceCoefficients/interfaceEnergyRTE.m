function out = interfaceEnergyRTE(matArray, modeIn, thetaDeg, N, opts)
%INTERFACEENERGYRTE  Run the RTE through one flat interface and log energy.
%
%   out = interfaceEnergyRTE(matArray, modeIn, thetaDeg, N, opts) launches a
%   finite plane-wave packet of N particles of type modeIn ('P', 'SV' or
%   'SH') at thetaDeg from the normal, inside matArray{1} (z < 0), lets
%   propagateParticleMultiMaterial carry it across the interface at z = 0
%   into matArray{2} (z > 0), and records the energy in each of the six
%   outgoing channels at every output time.
%
%   Unlike the direct MaterialInterface* calls the rest of this folder makes,
%   this actually runs the transport: the particles travel, hit the
%   interface, and are measured afterwards from where they ended up.
%
%   out.t          1 x nTime, output times (s)
%   out.tHitCenter time at which the packet centre reaches z = 0 (s)
%   out.tHitMax    time by which every particle has reached z = 0 (s)
%   out.E          nTime x 7, columns
%                    1 R_P   2 R_SV  3 R_SH   (reflected, still in medium 1)
%                    4 T_P   5 T_SV  6 T_SH   (transmitted into medium 2)
%                    7 total
%   out.final      1 x 7, energies averaged over the settled window
%   out.labels     channel names
%
%   opts (optional): .aperture, .zCenter, .nTime, .tEndFactor
%                    .frameFcn  handle P -> anything, called at .nFrames equally spaced steps;
%                               results returned in out.frames (cell) at times out.frameT
%
%   Called with no arguments it runs a self-check on a shale/sandstone pair
%   at 25 degrees and compares every channel against ZoeppritzSolid, which is
%   also the smallest useful example of how to call it:
%
%     m1  = makeInterfaceMaterial(3000,1500,2200);   % medium 1, z < 0
%     m2  = makeInterfaceMaterial(4500,2600,2600);   % medium 2, z > 0
%     out = interfaceEnergyRTE({m1,m2}, 'P', 25, 50000);
%
%   See also PROPAGATEPARTICLEMULTIMATERIAL, BOUNDARYREFLECTIONENERGYTESTS.

if nargin == 0
    out = selfCheck();
    return
end
if nargin < 5, opts = struct; end
def = struct('aperture', 20, 'zCenter', 80, 'nTime', 160, 'tEndFactor', 1.25, ...
             'frameFcn', [], 'nFrames', 40);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
end

mat1 = matArray{1};
nPlane = [0 0 1];

% Incident direction: travelling from z < 0 towards the interface at z = 0,
% matching the convention of the direct-call examples in this folder
% (dir = [sin(theta) 0 cos(theta)], normal n = [0 0 1]).
dIn = [sind(thetaDeg), 0, cosd(thetaDeg)];
dIn = dIn / norm(dIn);

% Finite aperture perpendicular to dIn, so the packet crosses over a spread
% of times and the energy history shows the transition rather than a step.
uv = null(dIn);
u = uv(:,1).';
v = uv(:,2).';
c1 = (rand(N,1) - 0.5) * opts.aperture;
c2 = (rand(N,1) - 0.5) * opts.aperture;
x = [0 0 -opts.zCenter] + c1 .* u + c2 .* v;

if any(x(:,3) >= 0)
    error(['Some particles start on or past the interface. ', ...
           'Increase opts.zCenter or reduce opts.aperture.']);
end

dir = repmat(dIn, N, 1);
[eSV, eSH] = svshBasisLocal(dir, nPlane);

switch upper(modeIn)
    case 'P',  isP = true(N,1);  perp = eSH;  vIn = mat1.vp;
    case 'SV', isP = false(N,1); perp = eSV;  vIn = mat1.vs;
    case 'SH', isP = false(N,1); perp = eSH;  vIn = mat1.vs;
    otherwise, error('Unknown incident mode "%s". Use P, SV or SH.', modeIn);
end

P = struct('d', 3, 'N', N, 'x', x, 'dir', dir, 'perp', perp, ...
           'p', isP, 'alive', true(N,1), 't', zeros(N,1), ...
           'matIdx', ones(N,1));

% Hit times: dir_z > 0 and z0 < 0, so z(t) = z0 + vIn*dir_z*t reaches 0 at
hitTimes   = -x(:,3) ./ (vIn * dir(:,3));
tHitCenter = opts.zCenter / (vIn * cosd(thetaDeg));
tHitMax    = max(hitTimes);

geometry = struct('dimension', 3, 'frame', 'cartesian');
geometry.bnd(1) = struct('dir', 3, 'val', 0, 'type', '', 'mat', [1 2]);

tVec = linspace(0, opts.tEndFactor * tHitMax, opts.nTime);
E = zeros(opts.nTime, 7);
E(1,:) = channelEnergies(P, nPlane);
withFrames = ~isempty(opts.frameFcn);
frameIdx = round(linspace(1, opts.nTime, opts.nFrames));
frames = cell(1, opts.nFrames);
if withFrames, for j = find(frameIdx == 1), frames{j} = opts.frameFcn(P); end, end
for it = 2:opts.nTime
    P = propagateParticleMultiMaterial(matArray, geometry, P, tVec(it));
    E(it,:) = channelEnergies(P, nPlane);
    if withFrames, for j = find(frameIdx == it), frames{j} = opts.frameFcn(P); end, end
end

settled = tVec > 1.10 * tHitMax;
out = struct('t', tVec, 'tHitCenter', tHitCenter, 'tHitMax', tHitMax, ...
             'E', E, 'final', mean(E(settled,:), 1), ...
             'labels', {{'R_P','R_SV','R_SH','T_P','T_SV','T_SH','total'}});
if withFrames, out.frames = frames; out.frameT = tVec(frameIdx); end
end

% -------------------------------------------------------------------------
function out = selfCheck()
% Smallest run that would fail if the transport, the channel bookkeeping or
% the polarization handling broke: one pair, one angle, every channel checked
% against ZoeppritzSolid and the total required to close at 1.
fprintf('interfaceEnergyRTE self-check: shale/sandstone, P at 25 deg\n');
N = 50000; theta = 25;
vp1 = 3000; vs1 = 1500; rho1 = 2200;
vp2 = 4500; vs2 = 2600; rho2 = 2600;
matArray = {makeInterfaceMaterial(vp1,vs1,rho1), makeInterfaceMaterial(vp2,vs2,rho2)};

out = interfaceEnergyRTE(matArray, 'P', theta, N);

j1 = linspace(0,90,361);
z  = MaterialClass.ZoeppritzSolid(j1, vp1,vs1,rho1, vp2,vs2,rho2);
expected = [interp1(j1,z.E_Rpp,theta),  interp1(j1,z.E_Rpsv,theta), 0, ...
            interp1(j1,z.E_Tpp,theta),  interp1(j1,z.E_Tpsv,theta), 0, 1];

tol = 5/sqrt(N);
fprintf('  %-6s %10s %10s %10s\n','channel','should','actual','|diff|');
for k = 1:7
    d = abs(expected(k) - out.final(k));
    fprintf('  %-6s %10.4f %10.4f %10.4f   %s\n', out.labels{k}, ...
        expected(k), out.final(k), d, passTag(d < tol));
end
assert(all(abs(expected - out.final) < tol), ...
    'interfaceEnergyRTE self-check failed against ZoeppritzSolid');
assert(abs(out.final(7) - 1) < 1e-12, ...
    'interfaceEnergyRTE lost or duplicated particles: total = %.6f', out.final(7));
fprintf('  self-check PASSED (tolerance %.4f at N = %d)\n', tol, N);
end

% -------------------------------------------------------------------------
function s = passTag(ok)
if ok, s = '[PASS]'; else, s = '[FAIL]'; end
end

% -------------------------------------------------------------------------
function E = channelEnergies(P, nPlane)
% Energy fractions in the six outgoing channels. P energy is a plain particle
% count; S energy splits between SV and SH by the squared projection of the
% polarization onto the SV/SH basis of the interface, exactly as
% BoundaryReflectionEnergyTests.m does it.
E = zeros(1,7);
side = {P.matIdx == 1, P.matIdx == 2};
for s = 1:2
    inP = side{s} &  P.p;
    inS = side{s} & ~P.p;
    E(3*(s-1)+1) = sum(inP) / P.N;
    if any(inS)
        [eSV, eSH] = svshBasisLocal(P.dir(inS,:), nPlane);
        aSV = dot(P.perp(inS,:), eSV, 2);
        aSH = dot(P.perp(inS,:), eSH, 2);
        E(3*(s-1)+2) = sum(aSV.^2) / P.N;
        E(3*(s-1)+3) = sum(aSH.^2) / P.N;
    end
end
E(7) = sum(E(1:6));
end

% -------------------------------------------------------------------------
function [eSV, eSH] = svshBasisLocal(d, nPlane)
nMat = repmat(nPlane(:).', size(d,1), 1);
eSH = cross(d, nMat, 2);
normSH = vecnorm(eSH, 2, 2);
% At exactly normal incidence the SV/SH split is arbitrary; pick a fixed
% in-plane pair rather than failing, since the two are then degenerate.
degenerate = normSH < 1e-12;
eSH(degenerate,:) = repmat([0 1 0], nnz(degenerate), 1);
normSH(degenerate) = 1;
eSH = eSH ./ normSH;
eSV = cross(eSH, d, 2);
eSV = eSV ./ vecnorm(eSV, 2, 2);
end
