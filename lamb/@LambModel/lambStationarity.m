function out = lambStationarity(M, z, opts)
%LAMBSTATIONARITY  Stationarity of equipartition under the local slab rule (gate G4.8).
%   out = lambStationarity(M, z)
%   out = lambStationarity(M, z, 'NBin', 12, 'NTheta', 40, 'NGLBin', 16)
%
%   Deterministic, no MC. For the uniform equipartition field
%   (energy per mode equal; body cell energy density u_a = w^2/(2 pi^2 c_a^3)/(4 pi) per solid
%   angle and polarisation; IW w/(2 pi c_I^2) per area), at each depth z (either medium) it
%   returns, per body cell (type b of the medium at z, sense in = toward z = 0 / out, theta
%   bin, phi-integrated), the energy removed R and injected J per unit volume and time.
%
%   Local rule (rule B; rule A, incoherent labels and fixed evanescent shares were rejected).
%   The diffuse field at z in medium m is the incoherent sum of all incident
%   states of both media. Each state at z = incident wave s_i (if from m) + evanescent waves
%   s_e + propagating outgoing waves c_o s_o (unit plane waves of m, scatteringState 'Away');
%   flux normalisation N_a = I_a cos(theta_a):
%     in-cell a  : [W(s_i) + f own_a + sum_o Re(conj c_oa X(s_o, s_i))] / N_a
%     out-cell b : W(s_b)/N_b + (1-f) own_b/N_b + sum_a Re(conj c_ba X(s_b, s_i,a + s_e,a))/N_a
%     tunnelling : evanescent waves transmitted from an incident state of the OTHER medium,
%                  [W(s_e) + sum_o Re(conj c_o X(s_o, s_e))]/N_a: owned by that state's in-cell,
%                  realised as an interface event (out.Rtun)
%   own_a = W(s_e) + 2 Re X(s_e, s_i) + sum_o Re(conj c_oa X(s_o, s_e)) (evanescent reflections),
%   f(k, z) such that the in-cell and out-cell removals of the type are equal (reciprocity
%   pairs emission-out with removal-in). X(y, v) = y' Q v, Q the total power form at z.
%
%   Emission (lambEmitLocal): channel ~ W_c, direction ~ c_n |A|^2 on the nodes x
%   72 azimuth bins, A = A_d (direct) + sum_t A_t (via the toward wave t). Coherent labels:
%   every cross term split half / half between its two waves; direct -> out cell b, via a
%   propagating t -> in cell t; evanescent-via terms of a same-medium channel b: share f to
%   out b, 1 - f to in b; of an other-medium channel: tunnelling emission, injected at the
%   interface into that out cell (out.Jbnd). Signed weights: out.neg = [negative, positive].
%
%   Checks: J = R per cell integrated over depth windows (the diffuse field has interference
%   fringes, so the balance is not pointwise in z); Jbnd = Rtun per type and bin integrated
%   over the other medium's depths; total balance per z (out.balance).
%
%   Rule (opts.Rule): 'flip' (default) = rule B' (evanescent share 1/2 everywhere + mirror flips
%   at rates fin (in -> out) / fout (out -> in) = max(-/+(A_in - A_out), 0); flips appear in R of
%   the source cell and in J of the mirrored cell); 'B' = the original rule B (auto share f).
%   Random media (het.law set): the same rule with Q_t replaced by the
%   nonlocal Gram of lambRandomGram (B(x, y) for tr(Q_t x y')), rates per unit volume (n = 1).
%
%   out.cells(m).types / .c / .edges{b}; out.R{m}{b}, out.J{m}{b} (2-by-bins-by-Nz, sense
%   1 = in, 2 = out; zero at depths of the other medium); out.Jbnd{m}{b}, out.Rtun{m}{b}
%   (bins-by-Nz) for the cells of medium m, from events at depths of the other medium.
arguments
    M LambModel
    z double
    opts.NBin   (1,1) double = 12
    opts.NTheta (1,1) double = 40
    opts.NGLBin (1,1) double = 16
    opts.NF     (1,1) double = 40            % f grid points per angular segment
    opts.Rule   char {mustBeMember(opts.Rule, {'B', 'flip'})} = 'flip'   % rule B' (default) or B
end
I = M.iface;  w = M.omega;
hasUp = ~isempty(I.upper);  nmed = 1 + hasUp;
cells = lambMediumCells(I);
for m = 1:nmed
    for b = 1:numel(cells(m).types)
        cells(m).edges{b} = uniquetol([linspace(0, pi/2, opts.NBin + 1), cells(m).thc{b}], 1e-9);
    end
end
nb = 72;  phb = 2*pi*((1:nb) - 0.5)/nb;
cI = M.iw(1).c;  EIW = w/(2*pi*cI^2);
Nz = numel(z);
for m = 1:nmed
    for b = 1:numel(cells(m).types)
        ne = numel(cells(m).edges{b}) - 1;
        out.R{m}{b} = zeros(2, ne, Nz);  out.J{m}{b} = zeros(2, ne, Nz);
        out.Jbnd{m}{b} = zeros(ne, Nz);  out.Rtun{m}{b} = zeros(ne, Nz);
    end
end
out.cells = cells;  out.z = z;  out.balance = zeros(Nz, 2);  out.neg = zeros(Nz, 2);
inT = {};  inC = [];  inM = [];  inB = [];
for m = 1:nmed
    inT = [inT, cells(m).types]; inC = [inC, cells(m).c]; %#ok<AGROW>
    inM = [inM, m*ones(1, numel(cells(m).types))];  inB = [inB, 1:numel(cells(m).types)]; %#ok<AGROW>
end

for iz = 1:Nz
    zz = z(iz);  m = 1 + (zz < 0);  o = 3 - m;  C = cells(m);
    assert(m == 1 || hasUp, 'lambStationarity: z < 0 needs an upper medium');
    psd = isfield(M.het{m}, 'law');  n0 = 1;  Qt = [];
    if ~psd, n0 = M.het{m}.n;  Qt = lambPowerForm(M, zz); end
    E = emissionTables(M, zz, cells);
    % --- f(k) of each type of m, on a grid clustered at the critical angles -------------
    F = cell(1, numel(C.types));
    for b = 1:numel(C.types)
        seg = unique([0, C.thc{b}, pi/2]);  thg = [];
        for q = 1:numel(seg) - 1
            ug = ((1:opts.NF) - 0.5)/opts.NF;  thg = [thg, seg(q) + (seg(q+1) - seg(q))*ug.^3]; %#ok<AGROW>
        end
        kg = (w/C.c(b))*sin(thg);  fg = 0.5*ones(size(kg));  fn = nan(size(kg));  fd = fn;
        if strcmp(opts.Rule, 'B')                                   % rule B' has f = 1/2
            for q = 1:numel(kg)
                Fq = formsAt(M, m, kg(q), zz, Qt, psd, opts.Rule);  fg(q) = Fq.f(b);
                if isfield(Fq, 'fnum'), fn(q) = Fq.fnum(b);  fd(q) = Fq.fden(b); end
            end
        end
        F{b} = struct('k', kg, 'f', fg, 'fnum', fn, 'fden', fd);
    end
    % --- removal: local cells, then tunnelling owned by the other medium's in-cells ----
    Rtot = 0;
    for b = 1:numel(C.types)
        eg = C.edges{b};
        for q = 1:numel(eg) - 1
            [th, wt] = sineGL(eg(q), eg(q+1), opts.NGLBin);
            for i = 1:numel(th)
                Fq = formsAt(M, m, (w/C.c(b))*sin(th(i)), zz, Qt, psd, opts.Rule);
                v = [Fq.rin(b); Fq.rout(b)];                       % scattering removal
                f = wt(i)*2*pi*sin(th(i))*cos(th(i))*uOf(w, C.c(b))*C.c(b)*n0;
                fl = [Fq.fin(b); Fq.fout(b)];                      % rule B' flips (0 for rule B)
                out.R{m}{b}(:, q, iz) = out.R{m}{b}(:, q, iz) + f*(v + fl);
                out.J{m}{b}(:, q, iz) = out.J{m}{b}(:, q, iz) + f*flipud(fl);   % in->out, out->in
                Rtot = Rtot + f*sum(v);
            end
        end
    end
    if hasUp
        Co = cells(o);
        for a = 1:numel(Co.types)
            eg = Co.edges{a};
            for q = 1:numel(eg) - 1
                [th, wt] = sineGL(eg(q), eg(q+1), opts.NGLBin);
                for i = 1:numel(th)
                    Fq = formsAt(M, m, (w/Co.c(a))*sin(th(i)), zz, Qt, psd, opts.Rule);
                    vt = Fq.rtun(a);
                    f = wt(i)*2*pi*sin(th(i))*cos(th(i))*uOf(w, Co.c(a))*Co.c(a)*n0;
                    out.Rtun{o}{a}(q, iz) = out.Rtun{o}{a}(q, iz) + f*vt;
                    Rtot = Rtot + f*vt;
                end
            end
        end
    end
    % --- injection: events of every incident state at z, plus the IW ------------------
    acc = initAcc(E);  JIW = 0;  neg = [0 0];
    for a = 1:numel(inT)
        ca = inC(a);
        brk = unique([0, asin(ca./inC(inC > ca*(1 + 1e-12))), pi/2]);
        for sgm = 1:numel(brk) - 1
            [th, wt] = sineGL(brk(sgm), brk(sgm+1), opts.NTheta);
            for i = 1:numel(th)
                if psd
                    Gr = lambRandomGram(M, M.het{m}, zz, w*sin(th(i))/ca);
                    St = Gr.states([Gr.states.med] == inM(a) & [Gr.states.type] == inB(a));
                    if isempty(St), continue, end                  % grazing within 1e-12: no field
                    rate = wt(i)*2*pi*sin(th(i))*uOf(w, ca)*ca*n0*bform(Gr.G, St.vec, St.vec)/St.I;
                    [acc, jiw, ng] = addEmissionPSD(Gr, St.vec, rate, acc, F, E.r);
                else
                    S = lamb.rte.scatteringState(I, inT{a}, sin(th(i))/ca, w, zz, 'Mode', 'state');
                    s = lamb.rte.equivalentSource(S, M.het{m}, w);
                    rate = wt(i)*2*pi*sin(th(i))*uOf(w, ca)*ca*n0*real(s'*Qt*s)/S.I;
                    [acc, jiw, ng] = addEmission(E, s, phb, rate, acc, F);
                end
                JIW = JIW + jiw;  neg = neg + ng;
            end
        end
    end
    if psd
        Gr = lambRandomGram(M, M.het{m}, zz, w/cI, 'IW', true);
        RIW = EIW*cI*n0*bform(Gr.G, Gr.iwVec, Gr.iwVec)/Gr.iwI;       % IW removal at z
        [acc, jiw, ng] = addEmissionPSD(Gr, Gr.iwVec, RIW, acc, F, E.r);
    else
        SI = lamb.rte.scatteringState(I, 'IW1', [], w, zz);
        sI = lamb.rte.equivalentSource(SI, M.het{m}, w);
        RIW = EIW*cI*n0*real(sI'*Qt*sI);                              % IW removal at z
        [acc, jiw, ng] = addEmission(E, sI, phb, RIW, acc, F);
    end
    JIW = JIW + jiw;  neg = neg + ng;
    % --- bin the flows in theta of the receiving cells ----------------------------------
    Jsum = JIW;
    for c = 1:numel(E.ch)
        ch = E.ch(c);
        Jsum = Jsum + sum(acc.d{c}) + sum(acc.v{c}(:)) + sum(acc.e{c});
        if ch.same                                                  % direct -> out cell b of m
            b = ch.cell;
            out.J{m}{b}(2, :, iz) = out.J{m}{b}(2, :, iz) + ...
                binFlows(ch.C.cells, acc.d{c}, (w/C.c(b))*sin(C.edges{b})/E.r).';
        else                                                        % tunnelling -> medium o
            b = ch.cell;
            out.Jbnd{o}{b}(:, iz) = out.Jbnd{o}{b}(:, iz) + ...
                binFlows(ch.C.cells, acc.e{c}, (w/cells(o).c(b))*sin(cells(o).edges{b})/E.r);
        end
        for t = 1:size(acc.v{c}, 2)                                 % via t -> in cell t of m
            out.J{m}{t}(1, :, iz) = out.J{m}{t}(1, :, iz) + ...
                binFlows(ch.C.cells, acc.v{c}(:, t), (w/C.c(t))*sin(C.edges{t})/E.r).';
        end
    end
    out.balance(iz, :) = [Jsum, Rtot + RIW];
    out.neg(iz, :) = neg;
end
end

% ==========================================================================================


function u = uOf(w, c)
u = w^2/(2*pi^2*c^3)/(4*pi);
end

function F = formsAt(M, m, k, z, Qt, psd, rule)
% rule-B removal values per cell at (k, z): .f(b) .rin(b) .rout(b) (types of m), .rtun(a)
% (types of the other medium); point scatterers via lambLocalForms, random media via the Gram
if ~psd
    L = lambLocalForms(M, m, k, z, Qt, rule);
    F.f = L.f;  F.rin = zeros(size(L.f));  F.rout = F.rin;  F.rtun = zeros(1, numel(L.Htun));
    F.fin = L.fin;  F.fout = L.fout;
    for b = find(~cellfun(@isempty, L.Hin)), F.rin(b) = real(trace(Qt*L.Hin{b}));  F.rout(b) = real(trace(Qt*L.Hout{b})); end
    for a = find(~cellfun(@isempty, L.Htun)), F.rtun(a) = real(trace(Qt*L.Htun{a})); end
    return
end
Gr = lambRandomGram(M, M.het{m}, z, k);
L = lambRandomForms(Gr, M, rule);
F.f = L.f;  F.rin = L.rin;  F.rout = L.rout;  F.rtun = L.rtun;  F.fnum = L.fnum;  F.fden = L.fden;
F.fin = L.fin;  F.fout = L.fout;
end

function fe = shareAt(Fb, k)
% evanescent share at k: ratio of the interpolated (smooth) numerator and denominator
n = interp1(Fb.k, Fb.fnum, k, 'linear', NaN);  d = interp1(Fb.k, Fb.fden, k, 'linear', NaN);
fe = n./d;  bad = ~isfinite(fe) | d == 0;
fe(bad) = interp1(Fb.k, Fb.f, k(bad), 'linear', 0.5);  fe(isnan(fe)) = 0.5;
end

function v = bform(G, x, y)
v = real(x.'*G*conj(y));
end

function [acc, jiw, neg] = addEmissionPSD(Gr, vec, rate, acc, F, r)
% addEmission for random media: phi-integrated piece-pair powers P(c, c', node) of the state
% vec, split by the same coherent labels
nu = Gr.nu;  X = cell(1, numel(Gr.ch));  tot = 0;
for c = 1:numel(Gr.ch)
    ch = Gr.ch(c);  npc = ch.npc;  N = size(ch.Q, 3);  nt = npc - 1;
    Vb = zeros(nu*npc, npc);
    for pc = 1:npc, Vb((1:nu) + nu*(pc - 1), pc) = vec; end
    P = pagemtimes(pagemtimes(Vb.', ch.Q), conj(Vb));             % npc x npc x N
    P = permute(P, [3 1 2]);                                       % N x npc x npc
    pr = ch.prop;  ev = ~pr;                                       % N x nt (via pieces)
    Pdd = real(P(:, 1, 1));  Pdt = real(P(:, 1, 2:end));  Pdt = reshape(Pdt, N, nt);
    Ptt = real(P(:, 2:end, 2:end));                                % N x nt x nt
    wd = Pdd + sum(Pdt.*pr, 2);
    we = sum(Ptt.*reshape(ev, N, nt, 1).*reshape(ev, N, 1, nt), [2 3]) + 2*sum(Pdt.*ev, 2) + ...
         sum(Ptt.*reshape(ev, N, nt, 1).*reshape(pr, N, 1, nt), [2 3]);
    wv = zeros(N, nt);
    for t = 1:nt
        oth = pr;  oth(:, t) = false;
        wv(:, t) = pr(:, t).*(Ptt(:, t, t) + Pdt(:, t) + sum(reshape(Ptt(:, :, t), N, nt).*ev, 2) + ...
                              sum(reshape(Ptt(:, t, :), N, nt).*oth, 2));
    end
    if ch.same                                                     % share f: out b / in b
        b = ch.cell;
        fe = shareAt(F{b}, ch.C.k*r);
        wd = wd + fe.*we;  wv(:, b) = wv(:, b) + (1 - fe).*we;  we = zeros(size(we));
    end
    X{c} = struct('d', wd, 'v', wv, 'e', we, 'all', [wd; wv(:); we]);
    tot = tot + sum(wd) + sum(wv(:)) + sum(we);
end
Xiw = 0;
for j = 1:numel(Gr.iwQ)
    na = size(Gr.iwQ{j}, 1);  Vb = zeros(na, 1);  npc = na/nu;
    for pc = 1:npc, Vb((1:nu) + nu*(pc - 1)) = vec; end
    Xiw = Xiw + real(Vb.'*Gr.iwQ{j}*conj(Vb));
end
tot = tot + Xiw;  g = rate/tot;  neg = [0 0];
if tot == 0, g = 0; end
for c = 1:numel(Gr.ch)
    acc.d{c} = acc.d{c} + g*X{c}.d;  acc.v{c} = acc.v{c} + g*X{c}.v;  acc.e{c} = acc.e{c} + g*X{c}.e;
    neg = neg + g*[sum(-min(X{c}.all, 0)), sum(max(X{c}.all, 0))];
end
jiw = g*Xiw;
end


function E = emissionTables(M, z, cells)
% the table set lambEmitLocal uses for a source at z, with the physical cell of each channel
I = M.iface;  w = M.omega;  m = 1 + (z < 0);
E.mz = ones(9, 1);  phys = @(nm) nm;
if m == 1
    T = M.T;  sd = T.solid;  E.h = (w/T.iface.lower.beta)*z;
elseif I.upper.isFluid
    T = M.T;  sd = T.fluid;  E.h = -(w/T.iface.lower.beta)*z;
else
    T = M.Tflip;  sd = T.solid;  E.h = -(w/T.iface.lower.beta)*z;  E.mz = [1 1 -1 1 1 1 1 -1 -1].';
    phys = @(nm) regexprep(regexprep(regexprep(nm, '1$', 'X'), '2$', '1'), 'X$', '2');
end
E.r = w/T.omega0;  E.sd = sd;
tm = cells(m).types;  to = {};  if numel(cells) == 2, to = cells(3 - m).types; end
for c = 1:numel(sd.body)
    nm = phys(sd.body(c).name);
    E.ch(c).C = sd.body(c);
    E.ch(c).same = any(strcmp(tm, nm));
    if E.ch(c).same, E.ch(c).cell = find(strcmp(tm, nm)); else, E.ch(c).cell = find(strcmp(to, nm)); end
end
end

function acc = initAcc(E)
for c = 1:numel(E.ch)
    N = numel(E.ch(c).C.k);  nt = size(E.ch(c).C.Tr, 2);
    acc.d{c} = zeros(N, 1);  acc.v{c} = zeros(N, nt);  acc.e{c} = zeros(N, 1);
end
end

function [acc, jiw, neg] = addEmission(E, s, phb, rate, acc, F)
% flows of one source s (physical frame, azimuth 0) with total emitted power 'rate'
nb = numel(phb);  r = E.r;
st = (s.*E.mz).*[r r r r^2 r^2 r^2 r^2 r^2 r^2].';                % table frame and units
A = zeros(9, nb);
for q = 1:nb, A(:, q) = lamb.rte.rotationMatrix(-phb(q))*st; end
X = cell(1, numel(E.ch));  tot = 0;
for c = 1:numel(E.ch)
    Ch = E.ch(c).C;  N = numel(Ch.k);  nt = size(Ch.Tr, 2);  h = E.h;
    Ad = (Ch.away.*exp(Ch.nuAway*h))*A;
    At = zeros(N, nb, nt);
    for t = 1:nt
        At(:, :, t) = (Ch.Tr(:, t).*reshape(Ch.tow(:, t, :), N, 9).*exp(-Ch.nuTow(:, t)*h))*A;
    end
    pr = reshape(abs(real(Ch.nuTow)) <= 1e-9*max(abs(Ch.nuTow), 1), N, 1, nt);
    Ae = sum(At.*~pr, 3);  Ap = At.*pr;  Sp = sum(Ap, 3);
    wd = abs(Ad).^2 + real(Ad.*conj(Sp));                          % direct owner
    wv = abs(Ap).^2 + real(Ad.*conj(Ap)) + real(Ae.*conj(Ap)) + real(Ap.*conj(Sp - Ap));
    we = abs(Ae).^2 + 2*real(Ad.*conj(Ae)) + real(Ae.*conj(Sp));   % evanescent owner
    if E.ch(c).same                                                % share f: out b / in b
        b = E.ch(c).cell;
        fe = interp1(F{b}.k, F{b}.f, Ch.k*r, 'linear', 0.5);
        wd = wd + fe.*we;  wv(:, :, b) = wv(:, :, b) + (1 - fe).*we;  we = zeros(size(we));
    end
    X{c} = struct('d', Ch.c.*sum(wd, 2), 'v', Ch.c.*reshape(sum(wv, 2), N, nt), ...
                  'e', Ch.c.*sum(we, 2), 'all', [wd(:); wv(:); we(:)].*repmat(Ch.c, nb*(2 + nt), 1));
    tot = tot + sum(X{c}.d) + sum(X{c}.v(:)) + sum(X{c}.e);
end
Xiw = 0;
for j = 1:numel(E.sd.iw)
    Pj = E.sd.iw(j);
    hj = Pj.fac*sum(Pj.q.*Pj.tow.*reshape(exp(-Pj.nuTow.'*E.h), 1, [], 1), 2);
    Xiw = Xiw + sum(abs(reshape(hj, 1, 9)*A).^2);
end
tot = tot + Xiw;  g = rate/tot;  neg = [0 0];
if tot == 0, g = 0; end                                            % no field here (e.g. SH into a fluid)
for c = 1:numel(E.ch)
    acc.d{c} = acc.d{c} + g*X{c}.d;  acc.v{c} = acc.v{c} + g*X{c}.v;  acc.e{c} = acc.e{c} + g*X{c}.e;
    neg = neg + g*[sum(-min(X{c}.all, 0)), sum(max(X{c}.all, 0))];
end
jiw = g*Xiw;
end




function v = binFlows(cells, f, kedges)
% distribute node flows f over k bins, each flow uniform in the cell's mapped variable t
nbin = numel(kedges) - 1;  v = zeros(nbin, 1);
for i = 1:size(cells, 1)
    if f(i) == 0, continue, end
    t1 = cells(i, 1);  t2 = cells(i, 2);  ka = cells(i, 3);  kb = cells(i, 4);
    tk = tOf(kedges, ka, kb, cells(i, 5));
    lo = max(min(tk(1:end-1), tk(2:end)), t1);  hi = min(max(tk(1:end-1), tk(2:end)), t2);
    v = v + f(i)*max(hi - lo, 0).'/(t2 - t1);
end
end

function t = tOf(k, ka, kb, typ)
k = min(max(k, min(ka, kb)), max(ka, kb));
if typ == 2
    t = (2*k - ka - kb)/(kb - ka);
else
    t = asin(sqrt(max(min((k.^2 - ka^2)/(kb^2 - ka^2), 1), 0)));
end
end

function [x, w] = sineGL(a, b, n)
% Gauss-Legendre on [a, b] after the sine map (sqrt kinks at both ends)
i = (1:n-1)';  bb = i./sqrt(4*i.^2 - 1);
[V, D] = eig(diag(bb, 1) + diag(bb, -1));
[t, k] = sort(diag(D));  wt = 2*V(1, k)'.^2;
x = a + (b - a)*(1 + sin(pi*t/2))/2;  w = (b - a)*pi/4*cos(pi*t/2).*wt;
end
