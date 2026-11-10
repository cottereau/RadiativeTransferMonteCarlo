function P = propagateParticleCylinder(M, P, dt)
%PROPAGATEPARTICLECYLINDER  Lamb-model particles in a solid cylinder with three free faces.
%   P = propagateParticleCylinder(M, P, dt)
%
%   M  lambRTEModel of a free surface (+ lambScatteringModel), with
%      M.geom = struct('R', R, 'h', h): solid 0 <= z <= h (z = depth below the top face), r <= R;
%      M.edge = lambEdgeTable(vp/vs) (90 degree edges r = R at z = 0 and z = h).
%   P  particles of initializeParticlesLamb (global frame = the top face's frame) plus P.face
%      (1 top, 2 bottom, 3 lateral: the face an interface-wave particle lives on).
%
%   Every face is the free surface of propagateParticleLamb in its own frame (depth along the
%   inward normal): body particles go to the frame of their nearest face (slab tables, overlap
%   at the edges by the nearest-face rule), interface-wave particles to their own face. Top and
%   bottom frames are global; the lateral one is the developed cylinder (s = R phi, -z, depth):
%   body particles that start the step there fly in the exact tangent-plane frame of their start
%   point (straight lines), rows emitted during the step and Rayleigh particles are wrapped back
%   (geodesics on the wall; curvature O(depth/R) ignored).
%   After the step: body particles that crossed a face other than their frame's are moved back to
%   it and crossed there (free-surface Zoeppritz, timing error <= dt); Rayleigh particles that
%   reached an edge are transmitted / reflected / radiated into P or S with M.edge at their
%   obliquity psi (radiation direction from the in-plane pattern, cosine along the edge
%   c sin(psi)/cR: k along the edge is conserved).
%   Before the step: a body particle whose in-plane line (normal to the edge) passes within
%   b/2 of an edge becomes a Rayleigh particle there, b = sum_faces sigma / sqrt(1 - a^2),
%   sigma = lambda_c pc(theta; psi), halved for S (only one S polarisation couples to the edge;
%   emitted and captured S are unpolarised here). Every such line meets a face first (the width
%   b lies outside the quarter space): it is judged where it meets the face (on the discs the
%   strip of hit points is an annulus: its width is stretched to the straight-edge area). Detailed balance
%   with the forward edge coefficients, checked per channel (dbcheck: P 0.9996, S 1.0000).
%   Cancellation (M.cancel) runs once per step on the face-local view, with face and position
%   cells (default Dxy = lambda_S, NPhi = 8).
%   Ceiling: dt must be short against the edge distances of interest (several crossings of
%   another face within one step are handled by up to 4 back-and-cross passes).
arguments
    M LambModel
    P struct
    dt (1,1) double {mustBePositive}
end
G = M.geom;
if ~isfield(P, 'face'), P.face = ones(P.N, 1); end
P.id = (1:P.N).';                                                  % row identity through the step
a = P.alive;  body = find(a & P.mode <= 2);  iw = a & P.mode >= 10;
fr = zeros(P.N, 1);  phf = atan2(P.x(:, 2), P.x(:, 1));
fr(body) = nearestFace(G, P.x(body, :));  fr(iw) = P.face(iw);
[fh, ph, th] = firstHit(M, P, body, dt);                           % frame of the face hit this step:
tHit = inf(P.N, 1);  tHit(body) = th;
cap = edgeCapture(M, P, dt, tHit);                                 % decided on the straight lines
if ~isempty(M.noIW) && M.noIW, cap.id = zeros(0, 1); end         % P/S reflections only: no edge capture
hit = fh > 0;  fr(body(hit)) = fh(hit);  phf(body(hit)) = ph(hit); % its plane passes the hit point
Mf = M;  Mf.cancel = [];  Mf.zMirror = [];
parts = {};
for f = 1:3
    idx = find(a & fr == f);
    if isempty(idx), continue, end
    [Q, fm] = toLocal(G, f, takeRows(P, idx), phf(idx));
    Q = propagateParticleLamb(Mf, Q, dt);
    parts{end+1} = toGlobal(G, f, Q, fm); %#ok<AGROW>
end
P = catRows(P, parts);
P = escapes(M, P);
[P, nIn] = applyCapture(M, P, cap, dt);
[P, nOut] = edgeIW(M, P);
if isfield(P, 'edgeLog'), P.edgeLog = P.edgeLog + [nOut nIn]; end   % energy Rayleigh -> body, body -> Rayleigh
if ~isempty(M.het) && ~isempty(M.cancel)
    P = cancelLocal(M, P);
end
end

% ------------------------------------------------------------------------------------------
function f = nearestFace(G, x)
r = hypot(x(:, 1), x(:, 2));
[~, f] = min([x(:, 3), G.h - x(:, 3), G.R - r], [], 2);
end

function [Q, fm] = toLocal(G, f, Q, phf)
% fm: lateral frames of the original rows (phi_f, z): tangent plane at azimuth phf (default:
% the particle's own); top / bottom frames are global
fm = [];
if nargin < 4, phf = atan2(Q.x(:, 2), Q.x(:, 1)); end
switch f
    case 2
        Q.x = [Q.x(:, 1), -Q.x(:, 2), G.h - Q.x(:, 3)];
        Q.dir = Q.dir.*[1 -1 -1];  Q.perp = Q.perp.*[1 -1 -1];
    case 3
        ph = phf(:);  fm = [ph, Q.x(:, 3)];
        [e1, e2, e3] = latFrame(ph);
        D = Q.x - [G.R*cos(ph), G.R*sin(ph), Q.x(:, 3)];
        Q.dir = [dot(Q.dir, e1, 2), dot(Q.dir, e2, 2), dot(Q.dir, e3, 2)];
        Q.perp = [dot(Q.perp, e1, 2), dot(Q.perp, e2, 2), dot(Q.perp, e3, 2)];
        Q.x = [G.R*ph + dot(D, e1, 2), -Q.x(:, 3), dot(D, e3, 2)];
end
end

function Q = toGlobal(G, f, Q, fm)
Q.face(:) = f;
if isfield(Q, 'rFace'), Q.rFace(Q.rFace < 0) = f; end              % reflected this step on face f
switch f
    case 2
        Q.x = [Q.x(:, 1), -Q.x(:, 2), G.h - Q.x(:, 3)];
        Q.dir = Q.dir.*[1 -1 -1];  Q.perp = Q.perp.*[1 -1 -1];
    case 3
        n0 = size(fm, 1);  n = Q.N;
        plane = false(n, 1);  plane(1:n0) = Q.mode(1:n0) <= 2;     % original body rows: exact frame
        ph = Q.x(:, 1)/G.R;  ph(plane) = fm(plane, 1);
        [e1, e2, e3] = latFrame(ph);
        X = [(G.R - Q.x(:, 3)).*cos(ph), (G.R - Q.x(:, 3)).*sin(ph), -Q.x(:, 2)];   % wrapped
        o = find(plane);
        if ~isempty(o)
            O = [G.R*cos(ph(o)), G.R*sin(ph(o)), fm(o, 2)];
            X(o, :) = O + (Q.x(o, 1) - G.R*ph(o)).*e1(o, :) + (Q.x(o, 2) + fm(o, 2)).*e2(o, :) + Q.x(o, 3).*e3(o, :);
        end
        Q.dir = Q.dir(:, 1).*e1 + Q.dir(:, 2).*e2 + Q.dir(:, 3).*e3;
        Q.perp = Q.perp(:, 1).*e1 + Q.perp(:, 2).*e2 + Q.perp(:, 3).*e3;
        Q.x = X;
end
end

function [e1, e2, e3] = latFrame(ph)
% lateral face at azimuth ph: e1 = theta hat, e2 = -z, e3 = inward normal -r hat
n = numel(ph);
e1 = [-sin(ph), cos(ph), zeros(n, 1)];  e2 = repmat([0 0 -1], n, 1);  e3 = [-cos(ph), -sin(ph), zeros(n, 1)];
end

% ------------------------------------------------------------------------------------------
function P = escapes(M, P)
% body particles beyond a face: back to the face along -dir, cross it there (no scattering)
G = M.geom;  Mflat = M;  Mflat.het = [];  Mflat.cancel = [];  Mflat.zMirror = [];
for pass = 1:4
    b = P.alive & P.mode <= 2;  r = hypot(P.x(:, 1), P.x(:, 2));
    ov = [-P.x(:, 3), P.x(:, 3) - G.h, r - G.R];                   % overshoot past each face
    [o, f] = max(ov, [], 2);
    e = find(b & o > 1e-9*G.R);
    if isempty(e), return, end
    c = chanSpeed(M, P.mode(e));
    for g = 1:3
        k = e(f(e) == g);  ck = c(f(e) == g);
        if isempty(k), continue, end
        d = P.dir(k, :);  x = P.x(k, :);
        switch g                                                   % path length back to the face
            case 1, s = -x(:, 3)./min(d(:, 3), -1e-12);
            case 2, s = (x(:, 3) - G.h)./max(d(:, 3), 1e-12);
            case 3
                bb = sum(x(:, 1:2).*d(:, 1:2), 2);  aa = max(sum(d(:, 1:2).^2, 2), 1e-24);
                cc = sum(x(:, 1:2).^2, 2) - G.R^2;
                s = (bb - sqrt(max(bb.^2 - aa.*cc, 0)))./aa;       % smaller root: the exit crossing
        end
        P.x(k, :) = x - s.*d;
        for q = 1:numel(k)                                         % ponytail: one call each (few per step)
            Q = takeRows(P, k([q q]));  Q.alive(2) = false;           % dead pad: the engine needs N > 1
            [Q, fm] = toLocal(G, g, Q);  Q.x(:, 3) = 0;
            Q = propagateParticleLamb(Mflat, Q, max(s(q)/ck(q), 1e-12));
            P = putRows(P, k(q), takeRows(toGlobal(G, g, Q, fm), 1));
        end
    end
end
end

function c = chanSpeed(M, mode)
c = M.iface.lower.alpha*ones(numel(mode), 1);  c(mode == 2) = M.iface.lower.beta;
end

% ------------------------------------------------------------------------------------------
function [P, n] = edgeIW(M, P)
% Rayleigh particles past an edge: back to it, then T / R / P / S with M.edge at psi
G = M.geom;  E = M.edge;  cR = M.iw(1).c;  n = 0;
iw = find(P.alive & P.mode == 11);
if isempty(iw), return, end
f = P.face(iw);  x = P.x(iw, :);  d = P.dir(iw, :);  r = hypot(x(:, 1), x(:, 2));
out = (f <= 2 & r > G.R) | (f == 3 & (x(:, 3) < 0 | x(:, 3) > G.h));
k = iw(out);
if isempty(k), return, end
for q = k.'
    x = P.x(q, :);  d = P.dir(q, :);  fA = P.face(q);
    if fA <= 2                                                     % back along the face to r = R
        bb = x(1:2)*d(1:2).';  aa = max(d(1:2)*d(1:2).', 1e-24);  cc = x(1:2)*x(1:2).' - G.R^2;
        s = (bb - sqrt(max(bb^2 - aa*cc, 0)))/aa;  x = x - s*d;
        top = fA == 1;
    else                                                           % back along the helix
        top = x(3) < 0;  ze = G.h*~top;
        s = (x(3) - ze)/d(3);  ph = atan2(x(2), x(1));
        [e1, ~, ~] = latFrame(ph);  dd1 = d*e1.';  ph = ph - s*dd1/G.R;
        [e1, ~, ~] = latFrame(ph);
        x = [G.R*cos(ph), G.R*sin(ph), ze];  d = dd1*e1 + [0 0 d(3)];
    end
    ph = atan2(x(2), x(1));  rh = [cos(ph), sin(ph), 0];  eh = [-sin(ph), cos(ph), 0];
    nz = [0 0 1 - 2*~top];                                         % inward normal of the top / bottom face
    if fA <= 2, nA = nz;  nB = -rh;  fB = 3;  else, nA = -rh;  nB = nz;  fB = 1 + ~top; end
    cp = max(-d*nB.', 0);  se = d*eh.';  psi = atan2d(abs(se), cp);
    co = [interp1(E.psi, E.T, psi), interp1(E.psi, E.R, psi), interp1(E.psi, E.P, psi), interp1(E.psi, E.S, psi)];
    ch = find(rand*sum(co) <= cumsum(co), 1);
    cpsi = cosd(psi);  spsi = sign(se)*sind(psi);
    P.x(q, :) = x;
    switch ch                                                      % then the rest of the path (s)
        case 1, P.dir(q, :) = cpsi*nA + spsi*eh;  P.face(q) = fB;  P.perp(q, :) = nB;  P = moveIW(G, P, q, s);
        case 2, P.dir(q, :) = cpsi*nB + spsi*eh;  P = moveIW(G, P, q, s);
        otherwise                                                  % body wave from the edge
            c = M.iface.lower.alpha;  pt = E.pP;  if ch == 4, c = M.iface.lower.beta;  pt = E.pS; end
            row = interp1(E.psi, pt, psi);  th = sampleBin(row, E.thEdges);
            al = min(c*abs(spsi)/cR, 1)*sign(spsi);
            dn = sqrt(1 - al^2)*(cos(th)*nB + sin(th)*nA) + al*eh;
            P.dir(q, :) = dn;  P.perp(q, :) = randPerp(dn);  P.mode(q) = ch - 2;
            P.x(q, :) = x + 1e-9*G.R*(nA + nB) + (c*s/cR)*dn;  n = n + P.E(q);   % converted energy
            if isfield(P, 'emLog'), P.emLog(ch - 2) = P.emLog(ch - 2) + P.E(q); end
    end
end
end

function th = sampleBin(row, ed)
c = cumsum(max(row, 0));  i = find(rand*c(end) <= c, 1);
th = ed(i) + rand*(ed(i + 1) - ed(i));
end

% ------------------------------------------------------------------------------------------
function cap = edgeCapture(M, P, dt, tHit)
% body particles whose in-plane line passes within b/2 of an edge (local, straight-edge
% geometry): their id, time, and the Rayleigh state at the edge (applied after the step by
% applyCapture). Inner half of b: the line reaches its closest approach t* <= dt before any
% face. Outer half: the line meets a face first; it is judged at that hit point (tHit <= dt),
% where it still aims at the edge, because after the reflection it no longer does.
G = M.geom;  E = M.edge;  cR = M.iw(1).c;
cap = struct('id', zeros(0, 1), 'ts', [], 'x', zeros(0, 3), 'dir', zeros(0, 3), 'perp', zeros(0, 3), 'face', [], 'mode', []);
dth = diff(E.thEdges(1:2));  bmax = 4*2*pi*M.iface.lower.alpha/M.omega;
b = find(P.alive & P.mode <= 2);
if isempty(b), return, end
c = chanSpeed(M, P.mode(b));  tH = tHit(b);  d = P.dir(b, :);
inner = edgeGeom(G, P.x(b, :), d, c);
outer = edgeGeom(G, P.x(b, :) + (c.*min(tH, dt)).*d, d, c);         % at the hit point (if any)
okI = inner.go & inner.ts >= 0 & inner.ts <= dt & inner.ts <= tH;
okO = outer.go & tH <= dt & outer.ts >= 0;
g = inner;
fn = fieldnames(g).';  for f = fn, g.(f{1})(okO, :) = outer.(f{1})(okO, :); end
g.ts(okO) = tH(okO);                                               % converted when it meets the face
ok = (okI | okO) & g.u < c*dt + bmax & g.v < c*dt + bmax;
sp = cR*abs(g.al)./c;  ok = ok & sp < 1;                           % sin(psi): k along the edge is conserved
if ~any(ok), return, end
lam = 2*pi*c/M.omega;  psi = asind(min(sp, 1));
isP = P.mode(b) == 1;
pc = zeros(numel(b), 2);
for q = find(ok).'
    pt = E.pP;  if ~isP(q), pt = E.pS; end
    row = interp1(E.psi, pt, psi(q));
    ib = @(t) min(max(floor(t/dth) + 1, 1), numel(row));
    pc(q, :) = lam(q)*[row(ib(g.th(q))), row(ib(pi/2 - g.th(q)))]/dth;      % sigma into top/bottom, lateral
    if ~isP(q), pc(q, :) = pc(q, :)/2; end                         % S: one of its two polarisations couples
end
bw = sum(pc, 2)./sqrt(max(1 - g.al.^2, realmin));
acc = g.del < bw/2;
fl = okO & g.v < 1e-9*G.R;                                         % hits on the top / bottom disc: the
ell = (bw/2).*g.u./max(g.del, realmin);                            % straight-edge strip u < ell becomes an
ellp = G.R - sqrt(max(G.R^2 - 2*G.R*ell, 0));                      % annulus; same area (per edge length)
acc(fl) = g.u(fl) < ellp(fl);                                      % up to ell' = R - sqrt(R^2 - 2 R ell)
cv = find(ok & acc);
for q = cv.'
    if okO(q), xe = P.x(b(q), :) + c(q)*tH(q)*d(q, :);  else, xe = P.x(b(q), :) + c(q)*g.ts(q)*d(q, :); end
    pe = atan2(xe(2), xe(1));  top = g.top(q);
    rq = [cos(pe), sin(pe), 0];  eq = [-sin(pe), cos(pe), 0];  nz = [0 0 2*top - 1];
    cps = sqrt(max(1 - sp(q)^2, 0));  sps = sign(g.al(q))*sp(q);
    if rand*sum(pc(q, :)) < pc(q, 1), fc = 1 + ~top;  dn = cps*(-rq) + sps*eq;  pn = nz;
    else,                             fc = 3;         dn = cps*nz + sps*eq;     pn = -rq; end
    cap.id(end+1, 1) = P.id(b(q));  cap.ts(end+1, 1) = g.ts(q);  cap.x(end+1, :) = [G.R*rq(1:2), G.h*~top];
    cap.dir(end+1, :) = dn;  cap.perp(end+1, :) = pn;  cap.face(end+1, 1) = fc;  cap.mode(end+1, 1) = P.mode(b(q));
end
end

function g = edgeGeom(G, x, d, c)
% local straight-edge geometry of lines (x, d, speed c) near the edge of their half: distances
% u (from the wall), v (from the flat face), in-plane velocity (du, dv), along-edge cosine al,
% closest-approach time ts, offset del, incoming angle th (from the flat face), go: toward both
r = hypot(x(:, 1), x(:, 2));  top = x(:, 3) < G.h/2;
u = G.R - r;  v = x(:, 3);  v(~top) = G.h - x(~top, 3);
ph = atan2(x(:, 2), x(:, 1));  rh = [cos(ph), sin(ph), zeros(numel(ph), 1)];  eh = [-sin(ph), cos(ph), zeros(numel(ph), 1)];
du = -c.*dot(d, rh, 2);  dv = c.*d(:, 3).*(2*top - 1);  al = dot(d, eh, 2);
w2 = max(du.^2 + dv.^2, realmin);  ts = -(u.*du + v.*dv)./w2;
g = struct('u', u, 'v', v, 'al', al, 'ts', ts, 'del', abs(u.*dv - v.*du)./sqrt(w2), ...
           'th', atan2(-dv, -du), 'top', top, 'go', du <= 0 & dv <= 0);
end

function [P, nIn] = applyCapture(M, P, cap, dt)
% the captured particles (by id; whatever happened to them this step) become Rayleigh particles
% at the edge, moved on for the rest of the step
nIn = 0;
if isempty(cap.id), return, end
[tf, loc] = ismember(cap.id, P.id);  tf = tf & loc > 0;
for j = find(tf).'
    q = loc(j);  if ~P.alive(q), continue, end
    nIn = nIn + P.E(q);
    if isfield(P, 'capLog'), P.capLog(cap.mode(j)) = P.capLog(cap.mode(j)) + P.E(q); end   % diagnostic: by channel
    P.mode(q) = 11;  P.face(q) = cap.face(j);  P.x(q, :) = cap.x(j, :);  P.dir(q, :) = cap.dir(j, :);  P.perp(q, :) = cap.perp(j, :);
    P = moveIW(M.geom, P, q, M.iw(1).c*(dt - cap.ts(j)));
end
end

function [fh, ph, tH] = firstHit(M, P, b, dt)
% face (0: none) first hit within dt by the straight flight of body rows b, its azimuth, and
% the hit time (inf beyond dt)
G = M.geom;  x = P.x(b, :);  d = P.dir(b, :);  c = chanSpeed(M, P.mode(b));  L = c*dt;
s = inf(numel(b), 3);
m = d(:, 3) < 0;  s(m, 1) = -x(m, 3)./d(m, 3);
m = d(:, 3) > 0;  s(m, 2) = (G.h - x(m, 3))./d(m, 3);
aa = sum(d(:, 1:2).^2, 2);  bb = sum(x(:, 1:2).*d(:, 1:2), 2);  cc = sum(x(:, 1:2).^2, 2) - G.R^2;
m = aa > 1e-24;  s(m, 3) = (-bb(m) + sqrt(max(bb(m).^2 - aa(m).*cc(m), 0)))./aa(m);
[sm, fh] = min(s, [], 2);  fh(sm > L) = 0;  tH = sm./c;  tH(sm > L) = inf;
xh = x + sm.*d;  ph = atan2(xh(:, 2), xh(:, 1));
end

function P = moveIW(G, P, q, len)
% move Rayleigh particle q by len along its face (straight on top / bottom, helix on the wall)
if P.face(q) <= 2, P.x(q, :) = P.x(q, :) + len*P.dir(q, :);  return, end
ph = atan2(P.x(q, 2), P.x(q, 1));  [e1, ~, ~] = latFrame(ph);
d1 = P.dir(q, :)*e1.';  d3 = P.dir(q, 3);
ph = ph + len*d1/G.R;  [e1, ~, ~] = latFrame(ph);
P.x(q, :) = [G.R*cos(ph), G.R*sin(ph), P.x(q, 3) + len*d3];  P.dir(q, :) = d1*e1 + [0 0 d3];
end

% ------------------------------------------------------------------------------------------
function P = cancelLocal(M, P)
% lambCancel on the face-local view (nearest face; Rayleigh: own face), face-separated cells
G = M.geom;  a = find(P.alive);
V = takeRows(P, a);  V.gid = a;
f = nearestFace(G, V.x);  iw = V.mode >= 10;  f(iw) = V.face(iw);
for g = 1:3
    k = find(f == g);
    if isempty(k), continue, end
    [Q, ~] = toLocal(G, g, takeRows(V, k));
    Q.x(:, 1) = Q.x(:, 1) + 1e6*g;                                 % separate the faces' cells
    V = putRows(V, k, Q);
end
o = M.cancel;  if ~isfield(o, 'Dxy'), o.Dxy = 2*pi/M.kb1; end
if ~isfield(o, 'NPhi'), o.NPhi = 8; end
if isfield(P, 'N0'), V.N0 = P.N0; end
args = namedargs2cell(o);
V = lambCancel(M, V, args{:});
E = V.E;
P = takeRows(P, V.gid);  P.E = E;  P.N0 = V.N0;
end

% ------------------------------------------------------------------------------------------
function Q = takeRows(P, rows)
Q = P;  n0 = P.N;
for f = fieldnames(P).'
    F = P.(f{1});
    if any(strcmp(f{1}, {'N', 'N0', 'edgeLog', 'capLog', 'emLog'})) || size(F, 1) ~= n0, continue, end
    Q.(f{1}) = F(rows, :);
end
Q.N = numel(rows);
end

function P = putRows(P, k, Q)
% rows 1..numel(k) of Q replace rows k of P; further rows of Q are appended
n0 = P.N;  nk = numel(k);
for f = fieldnames(P).'
    F = P.(f{1});
    if any(strcmp(f{1}, {'N', 'N0', 'edgeLog', 'capLog', 'emLog'})) || size(F, 1) ~= n0, continue, end
    if ~isfield(Q, f{1}), continue, end
    F(k, :) = Q.(f{1})(1:nk, :);
    P.(f{1}) = [F; Q.(f{1})(nk + 1:end, :)];
end
P.N = n0 + Q.N - nk;
end

function P = catRows(P, parts)
% the step's particles = the union of the face groups (dead rows dropped)
if isempty(parts), return, end
n0 = P.N;
for f = fieldnames(P).'
    F = P.(f{1});
    if any(strcmp(f{1}, {'N', 'N0', 'edgeLog', 'capLog', 'emLog'})) || size(F, 1) ~= n0, continue, end
    C = cellfun(@(Q) Q.(f{1})(Q.alive, :), parts, 'UniformOutput', false);
    P.(f{1}) = vertcat(C{:});
end
P.N = size(P.x, 1);
end

function p = randPerp(d)
p = cross(d, randn(size(d)), 2);
p = p./vecnorm(p, 2, 2);
end
