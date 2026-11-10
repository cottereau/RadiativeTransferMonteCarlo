function T = radiationTables(iface, opts)
%RADIATIONTABLES  Lamb radiation operator of an interface, for any point source and depth.
%   T = lamb.rte.radiationTables(iface)
%   T = lamb.rte.radiationTables(iface, 'NGauss', 64)
%
%   Everything radiated by a point source s (lamb.rte.sourceVector) at depth z is linear in s,
%   and its depth dependence is ANALYTIC, so no depth grid is stored. For each outgoing body
%   channel b and each quadrature node k_n,
%
%     D_b(k_n, phi=0; zeta) = sum_t Tr_bt(k_n) Pi^-_t(k_n) e^{-nu_t h}  +  Pi^+_b(k_n) e^{+nu_b h}
%
%   where t runs over the waves travelling toward the interface (P, SV, SH from a solid
%   source; P from a fluid source), Tr is the interface transfer for a unit incident wave
%   (Zoeppritz with phases), Pi^-/Pi^+ are the whole-space source amplitudes toward/away
%   from the interface (lamb.PointSource.directAmplitudes, per basis source), and h = |z|.
%   For each real pole j the residue is r (l.'S)/den with S linear in the incident
%   amplitudes, so the interface-wave amplitude at azimuth 0 is
%
%     h_j(zeta) = fac_j * sum_t q_jt Pi^-_t(k_p) e^{-nu_t(k_p) h},  dW_j/dphi = |h_j(phi).'*s|^2.
%
%   The tables are built at the reference frequency omega0 = beta_1 (k_beta1 = 1, so z = zeta).
%   The problem has no intrinsic length, so any omega follows by the scaling in
%   lamb.rte.channelQ / lamb.rte.channelAmplitudes.
%
%   Output T:
%     .iface .omega0 .bp .poles .leaky .peaks
%     .solid.body(c)  channels for a source in the lower solid (zeta >= 0)
%     .fluid.body(c)  channels for a source in the fluid (zeta < 0, solid-fluid only)
%       fields: .name ('P1' 'SV1' 'SH1' 'P2' 'SV2' 'SH2' 'Pf') .k .cells (quadrature
%               cells at omega0, lamb.quad.windowNodes) .c (node weight
%               w_n k_n E_m(k_n)/(4 pi^2)) .dOmega (per-solid-angle weight
%               k_m^2 cos(theta_n) E_m(k_n)/(4 pi^2)) .Tr (N-by-nt) .tow (N-by-nt-by-9) .nuTow (N-by-nt)
%               .away (N-by-9) .nuAway (N-by-1) .leakMask (N-by-numel(leaky))
%     .solid.iw(j), .fluid.iw(j)  real poles: .type .k .fac .q (1-by-nt) .tow (1-by-nt-by-9)
%               .nuTow (1-by-nt)
arguments
    iface (1,1) lamb.Interface
    opts.NGauss (1,1) double = 64
end
w0 = iface.lower.beta;
bp = iface.branchPoints(w0);
[poles, leaky] = lamb.findPoles(iface, w0);
T.iface = iface;  T.omega0 = w0;  T.bp = bp;  T.poles = poles;  T.leaky = leaky;
T.peaks = lamb.quad.peaksFor(iface, w0, poles, leaky);

T.solid = sideTables(iface, true, w0, bp, T.peaks, poles, leaky, opts.NGauss);
if strcmp(iface.kind, 'solid-fluid')
    T.fluid = sideTables(iface, false, w0, bp, T.peaks, poles, leaky, opts.NGauss);
end
end

% ------------------------------------------------------------------------------------------
function S = sideTables(iface, inSolid, w0, bp, peaks, poles, leaky, NG)
if inSolid
    med = iface.lower;  inc = {'P0', 'S0', 'C0'};  zs = +1;      % zs only fixes the side
else
    med = iface.upper;  inc = {'Pf0'};             zs = -1;
end
B = lamb.rte.basisSources(zs);
S.body = struct('name', {}, 'k', {}, 'cells', {}, 'c', {}, 'dOmega', {}, 'Tr', {}, 'tow', {}, 'nuTow', {}, ...
                'away', {}, 'nuAway', {}, 'leakMask', {});
for w = lamb.energy.outgoingWaves(iface, inSolid)
    km = w.med.wavenumber(w0, w.type);
    [k, wk, cells] = lamb.quad.windowNodes([0, bp(bp < km), km], NG, peaks);
    C.name = w.name;  C.k = k;  C.cells = cells;
    E    = lamb.energy.fluxDensity(w.med, k, w0, w.type);
    C.c  = wk.*k.*E/(4*pi^2);
    C.dOmega = km^2*sqrt(max(0, 1 - (k/km).^2)).*E/(4*pi^2);   % per unit solid angle
    C.Tr = transfer(iface, inc, k, w0, w.field);
    [C.tow, C.nuTow, awayAll, nuAll] = sourceAmps(B, med, k, w0, inSolid);
    if isempty(w.away)
        C.away = zeros(numel(k), 9);  C.nuAway = zeros(numel(k), 1);
    else
        idx = find(strcmp(w.away, {'Pd0', 'Sd0', 'Cd0'}));          % solid: P | SV | SH
        if isempty(idx), idx = 1; end                               % fluid: Uf0 (P)
        C.away = reshape(awayAll(:, idx, :), [], 9);  C.nuAway = nuAll(:, idx);
    end
    C.leakMask = false(numel(k), numel(leaky));
    for j = 1:numel(leaky)
        C.leakMask(:, j) = abs(k - leaky(j).a) < min(1e4*leaky(j).b, 1e-2*leaky(j).a);
    end
    S.body(end+1) = C;
end

S.iw = struct('type', {}, 'k', {}, 'fac', {}, 'q', {}, 'tow', {}, 'nuTow', {});
for p = poles
    kp = p.k;
    [r, lT, den] = lamb.poleNullVectors(iface, kp, w0);
    q = zeros(1, numel(inc));
    for t = 1:numel(inc)
        a = struct();
        for u = 1:numel(inc), a.(inc{u}) = double(u == t); end
        q(t) = lT*iface.rhs(kp, w0, a)/den;
    end
    names = iface.unknowns();  res = struct();
    for i = 1:numel(names), res.(names{i}) = r(i); end
    W2Dr = lamb.energy.interfaceFlux2D(iface, kp, w0, res);        % flux of the null vector
    [tow, nuTow] = sourceAmps(B, med, kp, w0, inSolid);
    S.iw(end+1) = struct('type', p.type, 'k', kp, 'fac', sqrt(kp/(2*pi)*W2Dr), ...
                         'q', q, 'tow', tow, 'nuTow', nuTow);
end
end

function Tr = transfer(iface, inc, k, w0, field)
% Scattered amplitude of one outgoing wave for unit incident amplitude of each incident wave.
N = numel(k);  Tr = zeros(N, numel(inc));
for t = 1:numel(inc)
    a = struct();
    for u = 1:numel(inc), a.(inc{u}) = double(u == t)*ones(N, 1); end
    x = iface.solve(a, k, w0);
    Tr(:, t) = x.(field);
end
end

function [tow, nuTow, away, nuAway] = sourceAmps(B, med, k, w0, inSolid)
% Whole-space amplitudes of the 9 basis sources at azimuth 0: toward the interface (tow) and
% away from it (away), with the vertical wavenumber of each wave type.
N = numel(k);
if inSolid
    epsTow = -1;                                   % interface above a solid source
    nu = [med.verticalWavenumber(k, w0, 'P'), med.verticalWavenumber(k, w0, 'S'), ...
          med.verticalWavenumber(k, w0, 'S')];
    nt = 3;
else
    epsTow = +1;                                   % interface below a fluid source
    nu = med.verticalWavenumber(k, w0, 'P');
    nt = 1;
end
tow = zeros(N, nt, 9);  away = zeros(N, nt, 9);
for i = 1:9
    [a1, a2, a3] = B{i}.directAmplitudes(med, k, 0, w0, epsTow);
    [b1, b2, b3] = B{i}.directAmplitudes(med, k, 0, w0, -epsTow);
    A = [a1, a2, a3];  Bw = [b1, b2, b3];
    tow(:, :, i) = A(:, 1:nt);  away(:, :, i) = Bw(:, 1:nt);
end
nuTow = nu;  nuAway = nu;
end
