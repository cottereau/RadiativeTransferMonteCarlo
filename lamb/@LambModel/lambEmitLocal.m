function P = lambEmitLocal(M, H, zs, xy, E, psd)
%LAMBEMITLOCAL  Emission of scattering events of the local slab transport.
%   P = lambEmitLocal(M, H, zs, xy, E)
%
%   H   9-by-9-by-Ne Hermitian source forms (physical frame, global azimuth; rank one s s' for
%       a plain source), e.g. the rule-B forms of lambLocalTables rotated to the particle's
%       azimuth: R H R' with R = lamb.rte.rotationMatrix(phi)
%   zs  1-by-Ne event depths;  xy  Ne-by-2 positions;  E  Ne-by-1 SIGNED energy to emit
%   P   particles (x, dir, perp, mode, side, t, E, alive, from); P.E is signed; P.from(i) is
%       the event that emitted row i: an event whose weights nearly cancel (sum|w| >> |sum w|)
%       emits n = ceil(sum|w|/|sum w|) (at most 64) independent particles, each with 1/n of the
%       energy factor, so that single weights do not grow
%
%   The emitted power pattern of H in the direction with amplitude row a is a H a' (Lamb
%   tables of the source side, as K1 (initializeParticlesLamb): nodes x 72 azimuth bins, the node's k
%   sampled in its cell, phi jittered in its bin). With A = A_d (direct) + sum_t A_t (via the
%   toward wave t) the pattern is split by the COHERENT labels of lambStationarity (validated,
%   gate G4.8): every cross term half / half between its two waves;
%     direct                     -> body particle heading away (out cell of that type)
%     via a propagating t        -> body particle of type t heading toward z = 0 (in cell)
%     evanescent-via, same medium -> share f(k, z) (loc(m).f) to the out cell, 1 - f to the in
%                                  cell of the channel's type
%     evanescent-via, other medium -> tunnelling: particle of that channel in the other medium,
%                                  at the interface, heading away
%     real interface-wave poles  -> IW particle;  leaky peaks of Pf -> leaky particle
%   Random medium: psd(e) per event (.Hc component form, .s 9 x nc x 3 component
%   sources at the event depth in the physical global frame, .K vertical wavenumbers, .k and .ph
%   incident horizontal wavenumber and azimuth, .W random-medium total of the form, .het): H is
%   then only the proposal; each drawn outcome o (channel, node, azimuth bin, label) is
%   reweighted by (w_psd(o)/W) / (proposal share of o), w_psd from the nonlocal pair weights
%   (lambRandomGram), so the emission follows the random-medium pattern in expectation.
%   The weights can be negative (1-3 % of the mass): one outcome is drawn per event with
%   probability |w|/sum|w| and carries energy E sum|w| sign(w) / sum(w), so the emitted energy
%   is right in expectation.
Ne = numel(zs);  zs = zs(:).';  E = E(:);
P = blank(0);  P.from = zeros(0, 1);
nc = 16;                                   % events per chunk: bounds the node x 72 x rank arrays
for a = 1:nc:Ne
    m = false(1, Ne);  m(a:min(a + nc - 1, Ne)) = true;
    for side = 1:2
        g = find(m & ((zs >= 0) == (side == 1)));
        if isempty(g), continue, end
        pg = [];  if nargin > 5 && ~isempty(psd), pg = psd(g); end
        Q = group(M, H(:, :, g), zs(g), xy(g, :), E(g), side, pg);
        Q.from = reshape(g(Q.from), [], 1);
        for f = fieldnames(P).', P.(f{1}) = [P.(f{1}); Q.(f{1})]; end
    end
end
P.N = numel(P.E);
end

% ------------------------------------------------------------------------------------------
function Q = group(M, H, zs, xy, E, side, psd)
I = M.iface;  w = M.omega;  Ne = numel(zs);
L = [];  if ~isempty(M.loc), L = M.loc(side); end                 % no slab: evanescent share 1/2
cells = lambMediumCells(I);  typesM = cells(side).types;  cM = cells(side).c;
typesO = {};  cO = [];  if numel(cells) == 2, typesO = cells(3 - side).types;  cO = cells(3 - side).c; end
D = eye(9);  phys = @(nm) nm;
if side == 1
    T = M.T;  sd = T.solid;  h = (w/T.iface.lower.beta)*zs;
elseif I.upper.isFluid
    T = M.T;  sd = T.fluid;  h = -(w/T.iface.lower.beta)*zs;
else
    T = M.Tflip;  sd = T.solid;  h = -(w/T.iface.lower.beta)*zs;  D = diag([1 1 -1 1 1 1 1 -1 -1]);
    phys = @(nm) regexprep(regexprep(regexprep(nm, '1$', 'X'), '2$', '1'), 'X$', '2');
end
r = w/T.omega0;  Sc = diag([r r r r^2 r^2 r^2 r^2 r^2 r^2]);
[R9, R72] = rotations();  nb = 72;  phb = 2*pi*((1:nb) - 0.5)/nb;
% --- eigencomponents of the forms, in table frame and units ------------------------------
V = zeros(9, 0);  lam = zeros(0, 1);  ev = zeros(0, 1);
for e = 1:Ne
    Ht = Sc*D*double(H(:, :, e))*D*Sc;  Ht = (Ht + Ht')/2;
    [U, S] = eig(Ht);  s = real(diag(S));  keep = abs(s) > 1e-8*max(abs(s));
    V = [V, U(:, keep)];  lam = [lam; s(keep)];  ev = [ev; e*ones(nnz(keep), 1)]; %#ok<AGROW>
end
Nq = numel(lam);
Sg = sparse(1:Nq, ev, lam, Nq, Ne);                                % component -> event, weight
A9 = permute(pagemtimes(R9, V), [1 3 2]);                          % 9 x 9(phi) x Nq
hq = h(ev);
% --- stage 1: azimuth-integrated label weights (9 azimuths: exact, degree <= 4 in phi) -----
out = struct('W', {}, 'ch', {}, 'lab', {}, 'kind', {}, 'j', {}, 'same', {});
F = cell(1, numel(sd.body));
for c = 1:numel(sd.body)
    C = sd.body(c);  N = numel(C.k);  nt = size(C.Tr, 2);
    nm = phys(C.name);  same = any(strcmp(typesM, nm));
    [wd, wv, we] = pieces(C, hq, A9);                              % N x Nq (x nt), phi-integrated
    agg = @(x) x*Sg;                                               % -> N x Ne
    wd = agg(wd);  we = agg(we);  wvE = zeros(N, Ne, nt);
    for t = 1:nt, wvE(:, :, t) = agg(wv(:, :, t)); end
    fe = [];
    if same                                                        % evanescent share f
        b = find(strcmp(typesM, nm));
        if isempty(L)
            fe = 0.5*ones(N, Ne);
        else
            za = min(max(abs(zs), abs(L.z(1))), abs(L.z(end)));    % clamp to the table depths
            kc = min(max(C.k*r, L.k(1)), L.k(end));                % clamp to the table's k range
            if isfield(L, 'fnum')                                  % random medium: f = num/den
                fe = interp2(abs(L.z), L.k, L.fnum(:, :, b), za, kc, 'linear')./ ...
                     interp2(abs(L.z), L.k, L.fden(:, :, b), za, kc, 'linear');
                bad = ~isfinite(fe);
                if any(bad(:)), f0 = interp2(abs(L.z), L.k, L.f(:, :, b), za, kc, 'linear');  fe(bad) = f0(bad); end
                fe(isnan(fe)) = 0.5;
            else
                fe = interp2(abs(L.z), L.k, L.f(:, :, b), za, kc, 'linear');   % N x Ne
            end
        end
        wd = wd + fe.*we;  wvE(:, :, b) = wvE(:, :, b) + (1 - fe).*we;  we = zeros(size(we));
    end
    lk = false(N, 1);  jl = zeros(N, 1);
    if strcmp(C.name, 'Pf') && ~isempty(C.leakMask)                 % leaky peaks: whole density
        for j = 1:size(C.leakMask, 2), lk(C.leakMask(:, j)) = true;  jl(C.leakMask(:, j)) = j; end
        out(end+1) = struct('W', C.c.*(wd + sum(wvE, 3) + we).*lk, 'ch', c, 'lab', 0, 'kind', 'leaky', 'j', jl, 'same', same); %#ok<AGROW>
        wd(lk, :) = 0;  wvE(lk, :, :) = 0;  we(lk, :) = 0;
    end
    if same, out(end+1) = struct('W', C.c.*wd, 'ch', c, 'lab', 0, 'kind', 'direct', 'j', [], 'same', same); end %#ok<AGROW>
    for t = 1:nt, out(end+1) = struct('W', C.c.*wvE(:, :, t), 'ch', c, 'lab', t, 'kind', 'via', 'j', [], 'same', same); end %#ok<AGROW>
    if ~same, out(end+1) = struct('W', C.c.*we, 'ch', c, 'lab', 0, 'kind', 'tunnel', 'j', [], 'same', same); end %#ok<AGROW>
    F{c} = fe;
end
for j = 1:numel(sd.iw)
    wi = sum(abs(pagemtimes(iwAmp(sd.iw(j), hq), A9)).^2, 2)*(2*pi/9);   % 1 x 1 x Nq
    out(end+1) = struct('W', reshape(wi, 1, Nq)*Sg, 'ch', 0, 'lab', 0, 'kind', 'iw', 'j', j, 'same', false); %#ok<AGROW>
end
Wall = cat(1, out.W);  off = [0, cumsum(arrayfun(@(o) size(o.W, 1), out))];
Tt = sum(Wall, 1);  Aa = sum(abs(Wall), 1);
cA = cumsum(abs(Wall), 1);
% --- stage 2: n_e independent outcomes per event (n_e = ceil(sum|w|/|sum w|), at most 64, so
% that no emitted weight grows much beyond |E|), azimuth from the outcome's 72-bin pattern ------
nE = min(max(ceil(Aa./abs(Tt) - 1e-9), 1), 64);
Q = blank(sum(nE));  Q.from = zeros(sum(nE), 1);  q = 0;
for e = 1:Ne
    for rep = 1:nE(e)
    q = q + 1;
    pk = find(cA(:, e) >= rand*cA(end, e), 1);
        ob = find(pk <= off(2:end), 1);  o = out(ob);  in = pk - off(ob);
        cq = find(ev == e);  Ae = permute(pagemtimes(R72, V(:, cq)), [1 3 2]);   % 9 x 72 x Re
        le = reshape(lam(cq), 1, 1, []);  he = h(e);
        if strcmp(o.kind, 'iw')
            w72 = sum(le.*abs(pagemtimes(iwAmp(sd.iw(o.j), he), Ae)).^2, 3);
        else
            C = sd.body(o.ch);  nt = size(C.Tr, 2);
            [wd, wv, we] = pieces(C, he, Ae, in);                      % 1 x 72 x Re (x nt), per bin
            wd = sum(le.*wd, 3);  we = sum(le.*we, 3);  wv = reshape(sum(le.*wv, 3), 1, nb, nt);
            if o.same
                b = find(strcmp(typesM, phys(C.name)));  fe = F{o.ch}(in, e);
                wd = wd + fe*we;  wv(:, :, b) = wv(:, :, b) + (1 - fe)*we;  we = 0*we;
            end
            switch o.kind
                case 'direct', w72 = wd;
                case 'via',    w72 = wv(:, :, o.lab);
                case 'tunnel', w72 = we;
                case 'leaky',  w72 = wd + sum(wv, 3) + we;
            end
        end
        ca = cumsum(abs(w72));  ib = find(ca >= rand*ca(end), 1);
        Q.E(q) = E(e)*Aa(e)/Tt(e)*sign(w72(ib))*ca(end)/abs(sum(w72))/nE(e);  Q.from(q) = e;
        if ~isempty(psd)                                           % random medium: reweight
            pe = Wall(pk, e)/Tt(e)*w72(ib)/sum(w72);               % proposal share of o
            fe = [];  if ~strcmp(o.kind, 'iw') && o.same, fe = F{o.ch}(in, e); end
            wo = psdOutcome(psd(e), o, in, ib, sd, he, T, D, Sc, w, r, fe, typesM, phys, side);
            Q.E(q) = Q.E(q)*(wo/psd(e).W)/pe;
        end
        ph = phb(ib) + (rand - 0.5)*(2*pi/nb);
        tv = [-sin(ph), cos(ph), 0];
        switch o.kind
            case 'iw'
                Q.x(q, :) = [xy(e, :), 0];  Q.dir(q, :) = [cos(ph), sin(ph), 0];  Q.perp(q, :) = [0 0 1];
                Q.mode(q) = 10 + o.j;  Q.side(q) = 1;  continue
            case 'leaky'
                Q.x(q, :) = [xy(e, :), 0];  Q.dir(q, :) = [cos(ph), sin(ph), 0];  Q.perp(q, :) = [0 0 1];
                Q.mode(q) = 20 + o.j(in);  Q.side(q) = 1;  continue
        end
        C = sd.body(o.ch);  kk = lamb.quad.cellSample(C.cells(in, :), rand)*r;
        nm = phys(C.name);
        switch o.kind
            case 'direct', ty = nm;  sd2 = side;  zz = zs(e);  toward = false;
            case 'via',    ty = typesM{o.lab};  sd2 = side;  zz = zs(e);  toward = true;
            case 'tunnel', ty = nm;  sd2 = 3 - side;  zz = (3 - 2*sd2)*1e-9*M.dslab;  toward = false;
        end
        if sd2 == side, cty = cM(strcmp(typesM, ty)); else, cty = cO(strcmp(typesO, ty)); end
        sT = min(kk*cty/w, 1);  cT = sqrt(1 - sT^2);
        away = 3 - 2*sd2;                                              % +1 lower side, -1 upper side
        dz = away*cT;  if toward, dz = -dz; end
        dr = [sT*cos(ph), sT*sin(ph), dz];
        if startsWith(ty, 'SV'), pp = cross(dr, tv);
        elseif startsWith(ty, 'SH'), pp = tv;
        else, pp = cross(dr, randn(1, 3));  pp = pp/norm(pp);
        end
        Q.x(q, :) = [xy(e, :), zz];  Q.dir(q, :) = dr;  Q.perp(q, :) = pp;
        Q.mode(q) = 1 + ~startsWith(ty, 'P');  Q.side(q) = sd2;
    end
end
end

function [wd, wv, we] = pieces(C, h, A, rows)
% coherent label pieces of channel C for components with depths h and rotated sources A
% (9 x nphi x Nq). All nodes, phi-integrated (N x Nq): pieces(C, h, A); one node, per phi
% bin (1 x nphi x Nq): pieces(C, h, A, row).
N = numel(C.k);  nt = size(C.Tr, 2);  nq = size(A, 3);  nphi = size(A, 2);
if nargin < 4, rows = 1:N; end
n = numel(rows);  h = reshape(h, 1, 1, []);
Ad = pagemtimes(C.away(rows, :).*exp(C.nuAway(rows).*h), A);       % n x nphi x nq
At = zeros(n, nphi, nq, nt);
for t = 1:nt
    dt = (C.Tr(rows, t).*reshape(C.tow(rows, t, :), n, 9)).*exp(-C.nuTow(rows, t).*h);
    At(:, :, :, t) = pagemtimes(dt, A);
end
pr = reshape(abs(real(C.nuTow(rows, :))) <= 1e-9*max(abs(C.nuTow(rows, :)), 1), n, 1, 1, nt);
Ae = sum(At.*~pr, 4);  Ap = At.*pr;  Sp = sum(Ap, 4);
wd = abs(Ad).^2 + real(Ad.*conj(Sp));
we = abs(Ae).^2 + 2*real(Ad.*conj(Ae)) + real(Ae.*conj(Sp));
wv = abs(Ap).^2 + real(Ad.*conj(Ap)) + real(Ae.*conj(Ap)) + real(Ap.*conj(Sp - Ap));
if nargin < 4                                                      % integrate over phi (exact)
    wd = reshape(sum(wd, 2), n, nq)*(2*pi/nphi);  we = reshape(sum(we, 2), n, nq)*(2*pi/nphi);
    wv = reshape(sum(wv, 2), n, nq, nt)*(2*pi/nphi);
end
end

function wo = psdOutcome(ps, o, in, ib, sd, he, T, D, Sc, w, r, fe, typesM, phys, side)
% random-medium weight of outcome o (node in, azimuth bin ib, label), in the units of the Gram:
% node weight x bin width x sum_{uv} Hc_uv sum_j T_uc conj(T_vc') (1/2)(g_uc + conj g_vc')
nb = 72;  phi = 2*pi*(ib - 0.5)/nb;
st = reshape(Sc*D*reshape(ps.s, 9, []), 9, size(ps.s, 2), []);  % table frame and units
Rb = lamb.rte.rotationMatrix(-phi);
y = abs(he)/(w/T.iface.lower.beta);                                % |z|
if strcmp(o.kind, 'iw')
    Pj = sd.iw(o.j);  nt = numel(Pj.nuTow);
    A = zeros(nt, 9);
    for t = 1:nt, A(t, :) = Pj.fac*Pj.q(t)*reshape(Pj.tow(1, t, :), 1, 9)*exp(-Pj.nuTow(t)*he); end
    e = -Pj.nuTow(:);  kap = Pj.k*r;  cw = 1;
else
    C = sd.body(o.ch);  nt = size(C.Tr, 2);
    A = zeros(1 + nt, 9);  A(1, :) = C.away(in, :)*exp(C.nuAway(in)*he);
    for t = 1:nt, A(1 + t, :) = C.Tr(in, t)*reshape(C.tow(in, t, :), 1, 9)*exp(-C.nuTow(in, t)*he); end
    e = [C.nuAway(in); -C.nuTow(in, :).'];  kap = C.k(in)*r;  cw = C.c(in);
end
Kh = ps.K(:).';  if side == 2, Kh = -Kh; end                       % half-line frame y = |z|
Kc = -1i*e*(w/T.iface.lower.beta);
dk = sqrt(max(ps.k^2 + kap^2 - 2*ps.k*kap*cos(phi - ps.ph), 0));
np = size(A, 1);  nc = size(st, 2);
Tm = reshape(A*Rb*reshape(st, 9, []), np, nc, []);               % piece x component x column
g = lambDepthCorr(ps.het, Kh + Kc, y, dk);                         % np x nc
P = zeros(np);
for c = 1:np
    for c2 = 1:np
        X = reshape(Tm(c, :, :), nc, []);  Y = reshape(Tm(c2, :, :), nc, []);
        P(c, c2) = sum(ps.Hc.*(X*Y').*((g(c, :).' + conj(g(c2, :)))/2), 'all');
    end
end
if strcmp(o.kind, 'iw')
    wo = real(sum(P, 'all'))*cw*(2*pi/nb);  return
end
pr = (abs(real(C.nuTow(in, :))) <= 1e-9*max(abs(C.nuTow(in, :)), 1));  ev = ~pr;
Pd = real(P(1, 2:end));  Pt = real(P(2:end, 2:end));
wd = real(P(1, 1)) + sum(Pd(pr));
we = sum(Pt(ev, ev), 'all') + 2*sum(Pd(ev)) + sum(Pt(ev, pr), 'all');
wv = zeros(1, nt);
for t = find(pr)
    oth = pr;  oth(t) = false;
    wv(t) = Pt(t, t) + Pd(t) + sum(Pt(ev, t)) + sum(Pt(t, oth));
end
if ~isempty(fe)
    b = find(strcmp(typesM, phys(C.name)));
    wd = wd + fe*we;  wv(b) = wv(b) + (1 - fe)*we;  we = 0;
end
switch o.kind
    case 'direct', v = wd;
    case 'via',    v = wv(o.lab);
    case 'tunnel', v = we;
    case 'leaky',  v = wd + sum(wv) + we;
end
wo = v*cw*(2*pi/nb);
end

function hj = iwAmp(Pj, h)
% interface-wave amplitude rows (1 x 9 x numel(h)) of a real pole
nt = numel(Pj.nuTow);  h = reshape(h, 1, []);
hj = Pj.fac*sum(Pj.q.*Pj.tow.*reshape(exp(-Pj.nuTow.'*h), 1, nt, 1, numel(h)), 2);
hj = reshape(hj, 1, 9, numel(h));
end

function [R9, R72] = rotations()
% rotationMatrix(-phi) at 9 exact-quadrature azimuths and at the 72 bin centres (9 x 9 x n)
persistent a b
if isempty(a)
    a = zeros(9, 9, 9);  for q = 1:9, a(:, :, q) = lamb.rte.rotationMatrix(-2*pi*(q - 1)/9); end
    b = zeros(9, 9, 72); for q = 1:72, b(:, :, q) = lamb.rte.rotationMatrix(-2*pi*(q - 0.5)/72); end
end
R9 = a;  R72 = b;
end


function P = blank(n)
P = struct('x', zeros(n, 3), 'dir', zeros(n, 3), 'perp', zeros(n, 3), 'mode', zeros(n, 1), ...
           'side', zeros(n, 1), 't', zeros(n, 1), 'E', zeros(n, 1), 'alive', true(n, 1));
end

function P = put(P, m, Q)
for f = fieldnames(Q).'
    if size(Q.(f{1}), 1) == nnz(m), P.(f{1})(m, :) = Q.(f{1}); end
end
end
