function M = lambScatteringModel(M, hetLower, hetUpper, opts)
%LAMBSCATTERINGMODEL  Point-scatterer heterogeneity and the scattering tables.
%   M = lambScatteringModel(M, hetLower, hetUpper)
%   M = lambScatteringModel(M, hetLower, hetUpper, 'dslab', d, 'zMirror', [zTop zBottom])
%
%   het*  struct: .n (scatterers per volume), .V, and ABSOLUTE contrasts .drho .dlam .dmu
%         (missing fields = 0). hetUpper is ignored for a free surface.
%   dslab slab half-thickness d (default M.dsrc = 3 lambda_S): inside |z| < d the local slab
%         transport (rule B) with the half-space scattering state (decision D1, confirmed by
%         G4.6); outside, whole-space Born (lambBulkScatter). 0 = no slab (whole space).
%   zMirror optional specular walls [zTop zBottom] (zTop < 0 < zBottom), e.g. for closed-box
%         equipartition tests; [] = open.
%
%   Adds to M:
%   .het{1,2} .dslab .zMirror
%   .bulk{side}   lambBulkScatter totals (.sigma, .pP, .rate)
%   .iwT(j)       interface-wave tables: .z .wz .nsig (n*sigma per unit length) .mu (rate per
%                 unit length) for the eigenstate, normalised to unit power per unit width
%   .loc          local slab tables (lambLocalTables); [] when dslab = 0
%   .cancel       lambCancel options for the signed Monte Carlo (struct(): defaults)
arguments
    M LambModel
    hetLower struct
    hetUpper struct = struct()
    opts.dslab   (1,1) double = NaN
    opts.zMirror double = []
    opts.LocOpts cell = {}                   % passed to lambLocalTables (e.g. {'NkSeg', 10, 'Nz', 32})
end
psd = isfield(hetLower, 'law');                                    % random medium
if psd
    if isempty(fieldnames(hetUpper)), hetUpper = hetLower; end
    M.het = {hetLower, hetUpper};
else
    M.het = {defaults(hetLower), defaults(hetUpper)};
end
M.dslab = opts.dslab;  if isnan(M.dslab), M.dslab = M.dsrc; end
M.zMirror = opts.zMirror;
I = M.iface;  w = M.omega;  hasUp = ~isempty(I.upper);

M.bulk{1} = lambBulkScatter('totals', I.lower, M.het{1}, w);
if hasUp, M.bulk{2} = lambBulkScatter('totals', I.upper, M.het{2}, w); end

% --- interface-wave tables --------------------------------------------------------------------
M.iwT = struct('z', {}, 'wz', {}, 'nsig', {}, 'mu', {}, 'psd', {});
for j = 1:numel(M.T.poles)
    kp = M.T.poles(j).k*w/M.T.omega0;
    L = I.lower;
    nuMin = min(real([L.verticalWavenumber(kp, w, 'P'), L.verticalWavenumber(kp, w, 'S')]));
    Zl = 20/nuMin;
    [x, wx] = gl(96);
    zz = Zl*(x + 1)/2;  ww = Zl*wx/2;  zz = zz.';  ww = ww.';
    if hasUp
        U = I.upper;  nus = real(U.verticalWavenumber(kp, w, 'P'));
        if ~U.isFluid, nus = min(nus, real(U.verticalWavenumber(kp, w, 'S'))); end
        Zu = 20/nus;  zu = -Zu*(x + 1)/2;  wu = Zu*wx/2;
        zz = [fliplr(zu.'), zz];  ww = [fliplr(wu.'), ww];
    end
    if psd, M.iwT(j) = iwRandom(M, j, zz, ww); continue, end
    S = lamb.rte.scatteringState(I, sprintf('IW%d', j), [], w, zz);
    s = zeros(9, numel(zz));
    for side = 1:1 + hasUp
        m = (zz >= 0) == (side == 1);
        Ss = S;  Ss.u = S.u(:, m);  Ss.eps = S.eps(:, :, m);
        s(:, m) = lamb.rte.equivalentSource(Ss, M.het{side}, w);
    end
    nsig = densityAt(M, zz).*powerAt(M, zz, s);                     % I = 1 (unit power per width)
    M.iwT(j) = struct('z', zz, 'wz', ww, 'nsig', nsig, 'mu', sum(ww.*nsig), 'psd', []);
end
M.loc = [];  if M.dslab > 0, M.loc = lambLocalTables(M, opts.LocOpts{:}); end
M.cancel = struct();                                               % lambCancel options (defaults)
end

function T = iwRandom(M, j, zz, ww)
% interface-wave tables for a random medium: nsig(z) = B(v, v) of the eigenstate v over the
% components (nonlocal Gram, n = 1; its depth integral is the exact half-space rate, G5.3), and
% per node the component data of its scattering events (lambEmitLocal psd)
kI = M.T.poles(j).k*M.omega/M.T.omega0;  nz = numel(zz);  nsig = zeros(1, nz);
X = struct('Hc', {}, 's', {}, 'K', {}, 'H9', {});
for iz = 1:nz
    Gr = lambRandomGram(M, M.het{1 + (zz(iz) < 0)}, zz(iz), kI, 'IW', true);
    v = Gr.iwVec;  Hc = v*v';  nsig(iz) = real(sum(Hc.*Gr.G, 'all'));
    H9 = zeros(9);  for c = 1:3, H9 = H9 + Gr.s(:, :, c)*Hc*Gr.s(:, :, c)'; end
    X(iz) = struct('Hc', Hc, 's', Gr.s, 'K', Gr.K, 'H9', H9);
end
T = struct('z', zz, 'wz', ww, 'nsig', nsig, 'mu', sum(ww.*nsig), 'psd', X);
end

function W = powerAt(M, z, s)
% total radiated power of sources s at depths z (lower solid, fluid, or mirrored upper solid)
W = zeros(1, numel(z));
lo = z >= 0;
if any(lo), W(lo) = lamb.rte.totalPower(M.T, M.kb1*z(lo), M.omega, s(:, lo)); end
up = ~lo;
if any(up)
    if strcmp(M.iface.kind, 'solid-solid')
        mz = [1 1 -1 1 1 1 1 -1 -1].';
        W(up) = lamb.rte.totalPower(M.Tflip, -(M.omega/M.iface.upper.beta)*z(up), M.omega, s(:, up).*mz);
    else
        W(up) = lamb.rte.totalPower(M.T, M.kb1*z(up), M.omega, s(:, up));
    end
end
end

function n = densityAt(M, z)
n = M.het{1}.n*(z >= 0) + M.het{2}.n*(z < 0);
end

function h = defaults(h)
for f = {'n', 'V', 'drho', 'dlam', 'dmu'}
    if ~isfield(h, f{1}), h.(f{1}) = 0; end
end
end

function [x, w] = gl(n)
i = (1:n-1)';  b = i./sqrt(4*i.^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, k] = sort(diag(D));  w = 2*V(1, k)'.^2;
end
