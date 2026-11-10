function P = propagateParticleLamb(M, P, dt)
%PROPAGATEPARTICLELAMB  Advance Lamb-model particles by dt.
%   P = propagateParticleLamb(M, P, dt)
%
%   Transport:
%   - body particles fly straight, event by event within the step; at z = 0 they are reflected /
%     transmitted with the Lamb real-k (Zoeppritz) energy coefficients;
%   - interface-wave particles (mode 10+j) move along z = 0; leaky ones (20+l) leak into the
%     fluid at rate 2 Im(k_L) per unit length (no scattering of leaky particles);
%   - optional specular walls M.zMirror = [zTop zBottom].
%   Scattering (only if lambScatteringModel was applied to M) — local slab transport, rule B,
%   validated by gates G4.4 (equipartition, 3 interfaces), G4.6 (Xu 2021) and G4.8:
%   - inside the slab |z| < d every body particle scatters where it is, at the local rate of
%     its cell (type, sense in / out, k, z): its own plane wave plus half of the in x out
%     interference terms of the diffuse field (tables M.loc, lambLocalTables). Events are drawn
%     along each flight segment from the cumulative |rate|; the event takes the sign of the
%     local rate (negative: the particle continues, a clone and a negative emission are added).
%   - emission from the cell's Hermitian form with the coherent signed labels (lambEmitLocal:
%     half-space Green's function, so interface waves appear here);
%   - tunnelling: an incoming particle reaching z = 0 may scatter in the evanescent field it
%     leaves in the other medium (events at depths of that medium) before crossing;
%   - outside the slab: whole-space Born scattering (lambBulkScatter), sampled along the
%     segments as well;
%   - K3: an interface-wave particle scatters at rate mu_j per unit length and is re-emitted
%     from its eigenstate at z* ~ n sigma(z), through lambEmitLocal;
%   - weights are signed: lambCancel (M.cancel) merges the signs within phase-space cells after
%     every step and keeps the particle count near its start.
%   Stability: at most about one scattering event per particle per step (the |E| spread grows
%   by ~10 % per event between cancellations); strong scattering needs a smaller dt.
%   M.dslab = 0: no slab, whole-space scattering everywhere (reference for identical media, G4.3).
%
%   See also LAMBRTEMODEL, LAMBSCATTERINGMODEL, LAMBLOCALTABLES, LAMBEMITLOCAL, LAMBCANCEL,
%   INITIALIZEPARTICLESLAMB, OBSERVELAMB.
%
%   References:
%     - Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Monte Carlo simulations of coupled
%       body- and Rayleigh-wave multiple scattering in elastic media, Geophys. J. Int. 228,
%       1213-1236. doi:10.1093/gji/ggab370
scat = ~isempty(M.het);
if scat
    P = transportLocal(M, P, dt);                                  % flight + slab / bulk scattering
else
    P = flight(M, P, dt);                                          % no heterogeneity
end
body = P.alive & P.mode <= 2;
if ~isempty(M.zMirror)
    P = mirrors(P, body, M.zMirror);
end
% --- interface waves ---------------------------------------------------------------------
for j = 1:numel(M.iw)
    m = P.alive & P.mode == 10 + j;
    P.x(m, :) = P.x(m, :) + (M.iw(j).c*dt).*P.dir(m, :);
    if scat && M.iwT(j).mu > 0
        ev = find(m & rand(P.N, 1) < 1 - exp(-M.iwT(j).mu*M.iw(j).c*dt));
        if ~isempty(ev), P = iwOperator(M, P, ev, j); end
    end
end
for j = 1:numel(M.leaky)
    Lk = M.leaky(j);
    m = find(P.alive & P.mode == 20 + j);
    te = -log(rand(numel(m), 1))/(Lk.rate*Lk.c);                    % leak time
    conv = m(te < dt);  stay = m(te >= dt);
    P.x(stay, :) = P.x(stay, :) + (Lk.c*dt).*P.dir(stay, :);
    if ~isempty(conv)
        tc = te(te < dt);
        P.x(conv, :) = P.x(conv, :) + (Lk.c*tc).*P.dir(conv, :);
        cf = sqrt(1 - Lk.sinf^2);
        P.dir(conv, :) = [Lk.sinf*P.dir(conv, 1:2), -cf*ones(numel(conv), 1)];   % up, into fluid
        P.perp(conv, :) = randPerp(P.dir(conv, :));
        P.mode(conv) = 1;  P.side(conv) = 2;
        P.x(conv, :) = P.x(conv, :) + (M.iface.upper.alpha*(dt - tc)).*P.dir(conv, :);
    end
end
if scat && ~isempty(M.cancel)              % signed MC: cancellation
    args = namedargs2cell(M.cancel);
    P = lambCancel(M, P, args{:});
end
P.t = P.t + dt;
end

% ------------------------------------------------------------------------------------------
function P = flight(M, P, dt)
% free flight with Zoeppritz crossings, event by event (no heterogeneity)
body = P.alive & P.mode <= 2;  rem = zeros(P.N, 1);  rem(body) = dt;  act = body;
for pass = 1:8
    if ~any(act), break, end
    c = speeds(M, P);  vz = c.*P.dir(:, 3);
    toward = act & ((P.side == 1 & vz < 0) | (P.side == 2 & vz > 0));
    te = inf(P.N, 1);  te(toward) = -P.x(toward, 3)./vz(toward);
    ev = act & te <= rem;  mv = act & ~ev;
    P.x(mv, :) = P.x(mv, :) + (c(mv, :).*rem(mv, :)).*P.dir(mv, :);  rem(mv) = 0;    % (x, :): 0x1 even when N = 1
    P.x(ev, :) = P.x(ev, :) + (c(ev, :).*te(ev, :)).*P.dir(ev, :);  rem(ev) = rem(ev) - te(ev);
    if any(ev), P.x(ev, 3) = 0;  P = crossInterface(M, P, find(ev)); end
    act = ev & rem > 0;
end
end

function P = iwOperator(M, P, idx, j)
% K3: an interface-wave particle scatters; re-emit from its eigenstate's equivalent source.
I = M.iface;  w = M.omega;  T = M.iwT(j);
n = numel(idx);  H = zeros(9, 9, n);  zst = zeros(1, n);  ps = [];
cdf = cumsum(T.wz.*abs(T.nsig));
for q = 1:n
    iz = find(cdf >= rand*cdf(end), 1, 'first');  zs = T.z(iz);
    if isfield(T, 'psd') && ~isempty(T.psd)                        % random medium
        d = P.dir(idx(q), :);  ph = atan2(d(2), d(1));  Rr = lamb.rte.rotationMatrix(ph);
        X = T.psd(iz);  sg = reshape(pagemtimes(Rr, reshape(X.s, 9, [])), size(X.s));
        H(:, :, q) = Rr*X.H9*Rr';  zst(q) = zs;
        ps = [ps, struct('Hc', X.Hc, 's', sg, 'K', X.K, 'k', w/M.iw(j).c, 'ph', ph, ...
                         'W', T.nsig(iz), 'het', M.het{1 + (zs < 0)})]; %#ok<AGROW>
        continue
    end
    S = lamb.rte.scatteringState(I, sprintf('IW%d', j), [], w, zs);
    d = P.dir(idx(q), :);
    s = lamb.rte.rotationMatrix(atan2(d(2), d(1)))*lamb.rte.equivalentSource(S, M.het{1 + (zs < 0)}, w);
    H(:, :, q) = s*s';  zst(q) = zs;
end
Q = lambEmitLocal(M, H, zst, P.x(idx, 1:2), P.E(idx), ps);
[P, ~] = applyEvents(P, zeros(P.N, 1), idx(:), ones(n, 1), Q);
end

function c = speeds(M, P)
c = zeros(P.N, 1);
body = P.alive & P.mode <= 2;
for ch = M.chan
    m = body & P.side == ch.side & P.mode == ch.mode;
    c(m) = ch.c;
end
end

function P = mirrors(P, body, zm)
m = body & P.side == 1 & P.x(:, 3) > zm(2);
P.x(m, 3) = 2*zm(2) - P.x(m, 3);  P.dir(m, 3) = -P.dir(m, 3);  P.perp(m, 3) = -P.perp(m, 3);
m = body & P.side == 2 & P.x(:, 3) < zm(1);
P.x(m, 3) = 2*zm(1) - P.x(m, 3);  P.dir(m, 3) = -P.dir(m, 3);  P.perp(m, 3) = -P.perp(m, 3);
end

function P = crossInterface(M, P, idx)
% Zoeppritz (Lamb real-k) split of the particles idx, which sit on z = 0.
d = P.dir(idx, :);
sinT = vecnorm(d(:, 1:2), 2, 2);
l = d(:, 1:2)./sinT;  l(sinT == 0, :) = repmat([1 0], nnz(sinT == 0), 1);
L = lambConventionMap(d, P.perp(idx, :), [0 0 1]);
isS = P.mode(idx) == 2;
typ = repmat({'P'}, numel(idx), 1);
typ(isS & rand(numel(idx), 1) < L.aSV.^2) = {'SV'};
typ(isS & strcmp(typ, 'P')) = {'SH'};
up = P.side(idx) == 2;                                % incident from the upper medium
key = cellfun(@(t, u) [t num2str(u)], typ, num2cell(up), 'UniformOutput', false);
for g = unique(key).'
    sel = find(strcmp(key, g{1}));
    ty = g{1}(1:end-1);  fromUp = g{1}(end) == '1';
    if fromUp && M.iface.upper.isFluid
        I = M.iface;  inc = 'Pf';  cInc = M.iface.upper.alpha;  map = @(n) sideOf(n, false);
    elseif fromUp
        I = M.ifaceFlip;  inc = [ty '1'];  cInc = speed(M.iface.upper, ty);  map = @(n) sideOf(n, true);
    else
        I = M.iface;  inc = [ty '1'];  cInc = speed(M.iface.lower, ty);  map = @(n) sideOf(n, false);
    end
    p  = sinT(sel)/cInc;
    if isfield(P, 'rMode')                            % last reflection (surface-interference weight)
        P.rMode(idx(sel)) = find(strcmp(ty, {'P', 'SV', 'SH'}));  P.rP(idx(sel)) = p;  P.rFace(idx(sel)) = -1;
    end
    rt = lamb.energy.rtCoefficients(I, inc, p, M.omega);
    names = setdiff(fieldnames(rt), {'p', 'sum'}, 'stable').';
    E = cell2mat(cellfun(@(n) rt.(n), names, 'UniformOutput', false));
    E(~isfinite(E)) = 0;  E = max(E, 0);
    cE = cumsum(E, 2);
    pick = sum(rand(numel(sel), 1).*cE(:, end) > cE, 2) + 1;
    for q = 1:numel(names)
        rows = sel(pick == q);
        if isempty(rows), continue, end
        n = names{q};  side = map(n);
        med = M.iface.lower;  if side == 2, med = M.iface.upper; end
        oty = regexprep(n, '[12f]$', '');
        cOut = speed(med, oty);
        sO = min(p(pick == q)*cOut, 1);  cO = sqrt(1 - sO.^2);
        ii = idx(rows);
        nd = [sO.*l(rows, :), (3 - 2*side)*cO];       % away from z = 0 into its side
        t  = [-l(rows, 2), l(rows, 1), zeros(numel(rows), 1)];
        switch oty
            case 'SV', pp = cross(nd, t, 2);  P.mode(ii) = 2;
            case 'SH', pp = t;                 P.mode(ii) = 2;
            otherwise, pp = randPerp(nd);      P.mode(ii) = 1;
        end
        P.dir(ii, :) = nd;  P.perp(ii, :) = pp;  P.side(ii) = side;
    end
end
end

function s = sideOf(name, flipped)
s = 1 + ~(name(end) == '1');                          % '1' -> side 1; '2' or 'f' -> side 2
if flipped, s = 3 - s; end
end

function c = speed(med, ty)
if strcmp(ty, 'P'), c = med.alpha; else, c = med.beta; end
end

function p = randPerp(d)
p = cross(d, randn(size(d)), 2);
p = p./vecnorm(p, 2, 2);
end

% ==========================================================================================
% Local slab transport (rule B; tables M.loc from lambLocalTables)
% ==========================================================================================
function P = transportLocal(M, P, dt)
% Body flight within the step, event by event: segments end at z = 0 or when the step is
% used; along each segment a local scattering event is drawn from the cumulative |rate| of
% the particle's cell (type, sense in / out, k) at its depths; the event takes the sign of
% the local rate. Positive: the particle is replaced by the emission of its cell's form.
% Negative: the particle continues, a clone is added and an emission of energy -E is made.
% At z = 0 an incoming particle first gets its tunnelling chance (events at depths of the
% other medium), then the Zoeppritz split. Emitted particles fly on with the rest of the step.
w = M.omega;
body = P.alive & P.mode <= 2;
rem = zeros(P.N, 1);  rem(body) = dt;  act = body;
for pass = 1:8                                  % ponytail: > 8 events per step: rest dropped
    if ~any(act), break, end
    i = find(act);  n = numel(i);
    c = speeds(M, P);  c = c(i);
    m = P.side(i);  dz = P.dir(i, 3);
    k = w*min(vecnorm(P.dir(i, 1:2), 2, 2), 1)./c;
    toward = (m == 1 & dz < 0) | (m == 2 & dz > 0);
    z0 = P.x(i, 3);  a0 = abs(z0);
    th = inf(n, 1);  th(toward) = -z0(toward)./(c(toward).*dz(toward));
    tseg = min(rem(i), th);  ds = c.*tseg;
    a1 = abs(z0 + c.*dz.*tseg);
    sense = 2 - toward;                                            % 1 in, 2 out
    [bq, wq] = components(M, P, i);                                % n x 2 table types, weights
    tau = zeros(n, 1);
    for q = 1:2 - 2*isempty(M.loc)                                 % no slab (dslab = 0): bulk only
        g = bq(:, q) > 0;
        if any(g), tau(g) = tau(g) + wq(g, q).*segDepth(M, m(g), bq(g, q), sense(g), k(g), a0(g), a1(g), ds(g)); end
    end
    d = M.dslab;                                                   % bulk part of the segment
    fo = zeros(n, 1);  lo = min(a0, a1);  hi = max(a0, a1);
    fo(lo >= d) = 1;
    cr = lo < d & hi > d;  fo(cr) = (hi(cr) - d)./(hi(cr) - lo(cr));
    rb = zeros(n, 1);
    for sd = 1:numel(M.bulk)
        for md = 1:2
            g = m == sd & P.mode(i) == md;  rb(g) = M.bulk{sd}.rate(md)./c(g);   % per unit path
        end
    end
    taub = rb.*fo.*ds;
    u = rand(n, 1);  ev = u < 1 - exp(-(tau + taub));
    isb = ev & rand(n, 1) < taub./max(tau + taub, realmin);        % the event is a bulk one
    if any(isb)                                                    % uniform on the bulk part
        j = find(isb);  xo = rand(numel(j), 1);
        s0 = zeros(numel(j), 1);                                   % path fraction where |z| >= d starts
        up = a1(j) >= a0(j);
        s0(up & lo(j) < d) = (d - a0(j(up & lo(j) < d)))./(a1(j(up & lo(j) < d)) - a0(j(up & lo(j) < d)));
        fr = zeros(numel(j), 1);
        fr(up) = s0(up) + xo(up).*(1 - s0(up));                    % moving away: bulk at the end
        dn = ~up;  fr(dn) = xo(dn).*fo(j(dn));                     % moving toward: bulk at the start
        P.x(i(j), :) = P.x(i(j), :) + (fr.*ds(j)).*P.dir(i(j), :);
        rem(i(j)) = rem(i(j)) - fr.*tseg(j);
        P = lambBulkScatter('scatter', M, P, i(j));
        ev(j) = false;                                             % handled; they fly on
        moved = false(n, 1);  moved(j) = true;
    else
        moved = false(n, 1);
    end
    % --- no event: to the end of the segment; interface arrivals ---------------------------
    ne = find(~ev & ~moved);  ne = ne(:);                          % column even for one particle
    P.x(i(ne), :) = P.x(i(ne), :) + ds(ne).*P.dir(i(ne), :);
    rem(i(ne)) = rem(i(ne)) - tseg(ne);
    h = ne(toward(ne) & th(ne) <= tseg(ne)*(1 + 1e-12));
    if ~isempty(h)
        P.x(i(h), 3) = 0;
        [P, rem, cross] = tunnelling(M, P, rem, i(h), k(h), bq(h, :), wq(h, :));
        if ~isempty(cross), P = crossInterface(M, P, cross); end
    end
    % --- scattering events -------------------------------------------------------------
    e = find(ev);  e = e(:);
    if ~isempty(e)
        ne2 = numel(e);  Hs = zeros(9, 9, ne2);  zs = zeros(1, ne2);  sg = zeros(ne2, 1);  ps = [];
        keep = true(ne2, 1);
        for q = 1:ne2
            j = e(q);  ts = -log(1 - rand*(1 - exp(-tau(j))));          % slab optical depth of the event
            [as, qc, sg(q), Hf, isFlip] = eventPoint(M, m(j), bq(j, :), wq(j, :), sense(j), k(j), a0(j), a1(j), ds(j), ts, tau(j));
            fr = 0;  if abs(a1(j) - a0(j)) > 1e-12, fr = (as - a0(j))/(a1(j) - a0(j)); else, fr = ts/tau(j); end
            fr = min(max(fr, 0), 1);
            P.x(i(j), :) = P.x(i(j), :) + (fr*ds(j))*P.dir(i(j), :);
            rem(i(j)) = rem(i(j)) - fr*tseg(j);
            if isFlip                                              % rule B' mirror flip: no emission
                P.dir(i(j), 3) = -P.dir(i(j), 3);  P.perp(i(j), 3) = -P.perp(i(j), 3);
                keep(q) = false;  continue
            end
            ph = atan2(P.dir(i(j), 2), P.dir(i(j), 1));  Rr = lamb.rte.rotationMatrix(ph);
            Hs(:, :, q) = Rr*Hf*Rr';  zs(q) = (3 - 2*m(j))*max(as, 1e-9*M.dslab);
            if isfield(M.loc(m(j)), 'psd') && ~isempty(M.loc(m(j)).psd)
                ps = [ps, psdEvent(M.loc(m(j)), sense(j), qb(qc, bq(j, :)), k(j), as, zs(q), ph)]; %#ok<AGROW>
            end
        end
        e = e(keep);  Hs = Hs(:, :, keep);  zs = zs(keep);  sg = sg(keep);
        if ~isempty(e)
            Q = lambEmitLocal(M, Hs, zs, P.x(i(e), 1:2), sg.*P.E(i(e)), ps);
            [P, rem] = applyEvents(P, rem, i(e), sg, Q);
        end
    end
    body = P.alive & P.mode <= 2;
    act = body & rem > 1e-12*dt;
end
end

function [bq, wq] = components(M, P, i)
% table type index and weight of the (up to two) components of each particle: P -> (1, 1);
% S -> SV and SH with aSV^2, aSH^2 (types 2, 3); a fluid has only P (type 1)
n = numel(i);  bq = zeros(n, 2);  wq = zeros(n, 2);
isP = P.mode(i) == 1;
bq(isP, 1) = 1;  wq(isP, 1) = 1;
s = find(~isP);
if ~isempty(s)
    L = lambConventionMap(P.dir(i(s), :), P.perp(i(s), :), [0 0 1]);
    bq(s, 1) = 2;  wq(s, 1) = L.aSV.^2;  bq(s, 2) = 3;  wq(s, 2) = L.aSH.^2;
end
end

function tau = segDepth(M, m, b, sense, k, a0, a1, ds)
% optical depth of |rate| along straight segments |z|: a0 -> a1 (path length ds), per particle
tau = zeros(numel(m), 1);
for g = unique([m, b, sense], 'rows').'
    sel = m == g(1) & b == g(2) & sense == g(3);
    L = M.loc(g(1));  C = cumAbs(L, g(3), g(2));  zc = (0:size(C, 2) - 1)*L.dz;
    kk = min(max(k(sel), L.k(1)), L.k(end));                       % clamp: grazing k beyond the last node
    ca = @(a) interp2(zc, L.k, C, min(a(sel), zc(end)), kk, 'linear');
    da = abs(a1(sel) - a0(sel));
    t = abs(ca(a1) - ca(a0))./max(da, eps).*ds(sel);
    flat = da < 1e-9*L.dz;                                         % grazing: rate x path
    if any(flat)
        as = a0(sel);  dss = ds(sel);
        rr = (abs(rateAt(L, g(3), g(2), kk, as)) + flipAt(L, g(3), g(2), kk, as)).*(as <= zc(end));   % no rate beyond the slab
        t(flat) = rr(flat).*dss(flat);
    end
    tau(sel) = t;
end
end

function [as, qc, sg, H, isFlip] = eventPoint(M, m, bq, wq, sense, k, a0, a1, ds, ts, tau)
% depth of the event (optical depth ts of |rate| + flip rate from a0 toward a1), component, rate
% sign, form; isFlip: the event is a rule-B' mirror flip (probability flip/(|rate| + flip))
L = M.loc(m);  zc = (0:numel(L.z))*L.dz;  k = min(max(k, L.k(1)), L.k(end));
if abs(a1 - a0) > 1e-9*L.dz
    G = zeros(1, numel(zc));
    for q = find(bq > 0), G = G + wq(q)*interp1(L.k, cumAbs(L, sense, bq(q)), k, 'linear'); end
    G = G + (0:numel(zc) - 1)*1e-15*max(G(end), 1);                % strictly increasing
    g0 = interp1(zc, G, min(a0, zc(end)));
    tgt = g0 + sign(a1 - a0)*ts*abs(a1 - a0)/ds;
    as = interp1(G, zc, min(max(tgt, G(1)), G(end)));
else
    as = a0;
end
r = zeros(1, 2);  rf = r;
for q = find(bq > 0), r(q) = wq(q)*rateAt(L, sense, bq(q), k, as);  rf(q) = wq(q)*flipAt(L, sense, bq(q), k, as); end
t = abs(r) + rf;
qc = find(rand*sum(t) <= cumsum(t), 1);  if isempty(qc), qc = find(bq > 0, 1); end
isFlip = rand*t(qc) < rf(qc);
sg = sign(r(qc));  if sg == 0, sg = 1; end
H = [];  if ~isFlip, H = formAt(L, sense, bq(qc), k, as); end
end

function b = qb(qc, bq)
b = bq(qc);
end

function ps = psdEvent(L, sense, b, k, a, zs, ph)
% random-medium data of one event (lambEmitLocal psd): component form of the cell (sense 1 in,
% 2 out, 3 tunnelling), Gram total, component sources at zs in the global frame (azimuth ph)
X = L.psd;  az = abs(L.z);  a = min(max(a, az(1)), az(end));  kc = min(max(k, L.k(1)), L.k(end));
ik = min(find(L.k <= kc, 1, 'last'), numel(L.k) - 1);  fk = (kc - L.k(ik))/(L.k(ik+1) - L.k(ik));
iz = min(find(az <= a, 1, 'last'), numel(az) - 1);     fz = (a - az(iz))/(az(iz+1) - az(iz));
switch sense, case 1, T = X.Hcin; case 2, T = X.Hcout; otherwise, T = X.Hctun; end
bl = @(A) double((1 - fk)*(1 - fz)*A(:, :, ik, iz) + fk*(1 - fz)*A(:, :, ik+1, iz) + ...
                 (1 - fk)*fz*A(:, :, ik, iz+1) + fk*fz*A(:, :, ik+1, iz+1));
Hc = bl(T(:, :, :, :, b));  G = bl(X.G);
if isfield(X, 'c')                                                 % exact K at k (sqrt kinks)
    nu = sqrt(complex(X.w^2./X.c(:).^2 - k^2));  Kk = [-X.sm*nu; X.sm*nu];
else
    Kk = (1 - fk)*X.K(:, ik) + fk*X.K(:, ik+1);
end                                                                % s ~ e^{i K z}: exact for the
s0 = (1 - fk)*X.sref(:, :, :, ik) + fk*X.sref(:, :, :, ik+1);      % unit plane-wave components
s = s0.*reshape(exp(1i*Kk*(zs - X.zref)), 1, [], 1);
Rr = lamb.rte.rotationMatrix(ph);
s = reshape(Rr*reshape(s, 9, []), size(s));
ps = struct('Hc', Hc, 's', s, 'K', Kk, 'k', k, 'ph', ph, 'W', real(sum(Hc.*G, 'all')), 'het', X.het);
end

function C = cumAbs(L, sense, b)
if sense == 1, C = L.CinA(:, :, b); else, C = L.CoutA(:, :, b); end
end

function r = rateAt(L, sense, b, k, a)
% signed local rate per unit path at (k, |z|), bilinear on the table nodes (clamped)
if sense == 1, R = L.rin(:, :, b); else, R = L.rout(:, :, b); end
az = abs(L.z);  k = min(max(k, L.k(1)), L.k(end));
r = interp2(az, L.k, R, min(max(a, az(1)), az(end)), k, 'linear');
r(isnan(r)) = 0;
end

function r = flipAt(L, sense, b, k, a)
% rule-B' mirror-flip rate per unit path at (k, |z|) (in -> out for sense 1), bilinear, clamped
r = 0;  if ~isfield(L, 'fin'), return, end
if sense == 1, R = L.fin(:, :, b); else, R = L.fout(:, :, b); end
az = abs(L.z);  k = min(max(k, L.k(1)), L.k(end));
r = interp2(az, L.k, R, min(max(a, az(1)), az(end)), k, 'linear');
r(isnan(r)) = 0;
end

function H = formAt(L, sense, b, k, a)
% form at (k, |z|), bilinear in the table nodes (clamped), double precision
if sense == 1, T = L.Hin; else, T = L.Hout; end
az = abs(L.z);  a = min(max(a, az(1)), az(end));  k = min(max(k, L.k(1)), L.k(end));
ik = min(find(L.k <= k, 1, 'last'), numel(L.k) - 1);  fk = (k - L.k(ik))/(L.k(ik+1) - L.k(ik));
iz = min(find(az <= a, 1, 'last'), numel(az) - 1);    fz = (a - az(iz))/(az(iz+1) - az(iz));
H = double((1 - fk)*(1 - fz)*T(:, :, ik, iz, b) + fk*(1 - fz)*T(:, :, ik+1, iz, b) + ...
           (1 - fk)*fz*T(:, :, ik, iz+1, b) + fk*fz*T(:, :, ik+1, iz+1, b));
end

function [P, rem, cross] = tunnelling(M, P, rem, idx, k, bq, wq)
% particles idx arriving at z = 0 (incoming): tunnelling scattering at depths of the other
% medium (events owned by the incoming cell), then the rest cross by Zoeppritz
cross = idx(:);
if numel(M.loc) < 2, return, end
n = numel(idx);  pt = zeros(n, 2);
for j = 1:n
    Lo = M.loc(3 - P.side(idx(j)));
    kj = min(max(k(j), Lo.k(1)), Lo.k(end));
    for q = find(bq(j, :) > 0), pt(j, q) = wq(j, q)*interp1(Lo.k, Lo.PtunA(:, bq(j, q)), kj, 'linear'); end
end
ev = find(rand(n, 1) < 1 - exp(-sum(pt, 2)));  ev = ev(:);
if isempty(ev), return, end
Hs = zeros(9, 9, numel(ev));  zs = zeros(1, numel(ev));  sg = zeros(numel(ev), 1);  ps = [];
for q = 1:numel(ev)
    j = ev(q);  so = 3 - P.side(idx(j));  Lo = M.loc(so);
    qc = find(rand*sum(pt(j, :)) <= cumsum(pt(j, :)), 1);  b = bq(j, qc);
    kj = min(max(k(j), Lo.k(1)), Lo.k(end));
    Crow = interp1(Lo.k, Lo.CtunA(:, :, b), kj, 'linear');  zc = (0:numel(Lo.z))*Lo.dz;
    Crow = Crow + (0:numel(zc) - 1)*1e-15*max(Crow(end), 1);
    as = interp1(Crow, zc, rand*Crow(end));
    az = abs(Lo.z);  a = min(max(as, az(1)), az(end));
    qv = interp2(az, Lo.k, Lo.qtun(:, :, b), a, kj, 'linear');  sg(q) = sign(qv);  if sg(q) == 0, sg(q) = 1; end
    T = Lo.Htun;  kk = min(max(k(j), Lo.k(1)), Lo.k(end));
    ik = min(find(Lo.k <= kk, 1, 'last'), numel(Lo.k) - 1);  fk = (kk - Lo.k(ik))/(Lo.k(ik+1) - Lo.k(ik));
    iz = min(find(az <= a, 1, 'last'), numel(az) - 1);        fz = (a - az(iz))/(az(iz+1) - az(iz));
    Hf = double((1 - fk)*(1 - fz)*T(:, :, ik, iz, b) + fk*(1 - fz)*T(:, :, ik+1, iz, b) + ...
                (1 - fk)*fz*T(:, :, ik, iz+1, b) + fk*fz*T(:, :, ik+1, iz+1, b));
    ph = atan2(P.dir(idx(j), 2), P.dir(idx(j), 1));  Rr = lamb.rte.rotationMatrix(ph);
    Hs(:, :, q) = Rr*Hf*Rr';  zs(q) = (3 - 2*so)*max(as, 1e-9*M.dslab);
    if isfield(Lo, 'psd') && ~isempty(Lo.psd), ps = [ps, psdEvent(Lo, 3, b, k(j), as, zs(q), ph)]; end %#ok<AGROW>
end
Q = lambEmitLocal(M, Hs, zs, P.x(idx(ev), 1:2), sg.*P.E(idx(ev)), ps);
pos = sg > 0;
cross = setdiff(idx(:), idx(ev(pos)));                             % removed ones do not cross
[P, rem] = applyEvents(P, rem, idx(ev), sg, Q);
cross = cross(P.alive(cross) & P.mode(cross) <= 2 & P.x(cross, 3) == 0);
end

function [P, rem] = applyEvents(P, rem, k, sg, Q)
% Events k (signs sg) and their emission Q (lambEmitLocal; Q.from = event of each row, several
% rows per event when its weights nearly cancel). Positive event: the particle is replaced by
% its first emitted row, the other rows are appended. Negative event: the particle stays, a
% clone and all its emitted rows are appended. New rows keep the event particle's remaining time.
k = k(:);  sg = sg(:);  from = Q.from(:);
[~, first] = unique(from, 'first');                                % first row of each event
ev1 = from(first);
isFirst = false(numel(from), 1);  isFirst(first) = true;
pos = sg(ev1) > 0;
P = setRows(P, k(ev1(pos)), Q, first(pos));                        % replacements
extra = find(~(isFirst & sg(from) > 0));                          % every other emitted row
clones = k(sg < 0);
P = appendRows(P, P, clones);
P = appendRows(P, Q, extra);
rem = [rem; rem(clones); rem(k(from(extra)))];
end

function P = setRows(P, k, Q, q)
for f = {'x', 'dir', 'perp', 'mode', 'side', 'E'}
    P.(f{1})(k, :) = Q.(f{1})(q, :);
end
P.alive(k) = true;
end

function P = appendRows(P, Q, q)
% append rows q of Q to P (fields of P with one row per particle; missing ones zero)
n0 = P.N;  nq = numel(q);
if nq == 0, return, end
for f = fieldnames(P).'
    F = P.(f{1});
    if size(F, 1) ~= n0 || strcmp(f{1}, 'N'), continue, end
    if isfield(Q, f{1}) && size(Q.(f{1}), 1) >= max(q), add = Q.(f{1})(q, :);
    else, add = zeros(nq, size(F, 2), 'like', F); end
    P.(f{1}) = [F; add];
end
P.N = n0 + nq;
end
