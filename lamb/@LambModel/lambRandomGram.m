function Gr = lambRandomGram(M, het, z, k, opts)
%LAMBRANDOMGRAM  Random-medium power forms of the near-interface components at (k, z).
%   Gr = lambRandomGram(M, het, z, k)
%   Gr = lambRandomGram(M, het, z, k, 'IW', true, 'Nphi', 64)
%
%   The random-medium replacement of the 9-by-9 total power form Q_t of rule B (lambLocalForms,
%   lambStationarity), for scatterers at depth z (medium m of z; het its random medium, with
%   .law, see lambDepthCorr) and horizontal wavenumber k (azimuth 0).
%   Components u (columns of every vector below): the incident unit wave of each type of m
%   propagating at k, then the outgoing unit wave of each type of m (evanescent included);
%   s_u ~ e^{i K_u z}. Every scattered term T_{u,c} = A_c R(-phi) s_u (outgoing piece c:
%   direct, via toward wave t) is paired with the nonlocal weight (1/2)[g(K_alpha) + conj
%   g(K_beta)] (lambDepthCorr; exact in the depth integral, rule B as a -> 0).
%     Gr.G         nu-by-nu Hermitian Gram: B(x, y) = x.' G conj(y) is the scattered power
%                  (all channels and interface waves) of the component combinations x, y; for a
%                  point scatterer B(x, y) = s_y' Q_t s_x (tr(Q_t x y') of lambLocalForms)
%     Gr.ch(c)     body channels: .Q (na-by-na-by-N, alpha = u + nu (piece - 1), phi-integrated,
%                  node weight C.c included), .npc, .prop (N-by-nt, via piece propagating),
%                  .C (table channel), .name (physical), .same (channel of medium m), .cell
%     Gr.iwQ{j}    interface-wave channels (na-by-na)
%     Gr.states(i) every incident body state (both media) propagating at k: .med .type .vec
%                  (state = vec over the components), .si .se (incident / evanescent parts),
%                  .c (outgoing-wave coefficients), .N = I cos(theta), .I
%     Gr.iwVec     with 'IW', true: the interface-wave eigenstate (k = its pole) over the
%                  components, and Gr.iwI its normalisation (power per unit width)
%     Gr.types, Gr.prop (types of m propagating at k), Gr.iInc, Gr.iOut (component indices)
%     Gr.s (9 x nu x 3 component sources at z, physical frame), Gr.K (s ~ e^{i K z})
arguments
    M LambModel
    het struct
    z (1,1) double
    k (1,1) double
    opts.IW (1,1) logical = false
    opts.Nphi (1,1) double = 64
end
I = M.iface;  w = M.omega;  p = k/w;
cells = lambMediumCells(I);  m = 1 + (z < 0);  sm = 3 - 2*m;  Cm = cells(m);  nq = numel(Cm.types);
med = I.lower;  if m == 2, med = I.upper; end
L = covColumns(med, het);
fld = @(X) [X.u(:); X.eps(:)];
% --- components: sources (9 x nu x 3 covariance columns), K, fields -------------------------
Gr.types = Cm.types;  Gr.prop = k < w./Cm.c*(1 - 1e-12);
nuv = sqrt(complex(w^2./Cm.c.^2 - k^2));
St = {};  K = [];
for b = find(Gr.prop)
    St{end+1} = lamb.rte.scatteringState(I, Cm.types{b}, p, w, z, 'Mode', 'incident'); %#ok<AGROW>
    K(end+1) = -sm*nuv(b); %#ok<AGROW>
end
Gr.iInc = zeros(1, nq);  Gr.iInc(Gr.prop) = 1:nnz(Gr.prop);
U = zeros(12, nq);
for q = 1:nq
    St{end+1} = lamb.rte.scatteringState(I, Cm.types{q}, p, w, z, 'Mode', 'incident', 'Away', true); %#ok<AGROW>
    K(end+1) = sm*nuv(q); %#ok<AGROW>
    U(:, q) = fld(St{end});
end
Gr.iOut = nnz(Gr.prop) + (1:nq);  nu = numel(St);
s = zeros(9, nu, 3);
for u = 1:nu
    for j = 1:3
        hj = struct('V', 1, 'n', 1, 'dlam', L(1, j), 'dmu', L(2, j), 'drho', L(3, j));
        s(:, u, j) = lamb.rte.equivalentSource(St{u}, hj, w);
    end
end
Gr.s = s;  Gr.K = K;                                               % physical frame, s ~ e^{i K z}
% --- states over the components -------------------------------------------------------------
un = vecnorm(U);  solve = @(t) ((U./un)\t)./un.';
Gr.states = struct('med', {}, 'type', {}, 'vec', {}, 'si', {}, 'se', {}, 'c', {}, 'N', {}, 'I', {});
for mm = 1:numel(cells)
    for a = 1:numel(cells(mm).types)
        ca = cells(mm).c(a);
        if k >= w/ca*(1 - 1e-12), continue, end
        Ss = lamb.rte.scatteringState(I, cells(mm).types{a}, p, w, z, 'Mode', 'state');
        tgt = fld(Ss);  si = zeros(nu, 1);
        if mm == m, si(Gr.iInc(a)) = 1;  tgt = tgt - fld(St{Gr.iInc(a)}); end
        c = solve(tgt);
        c(abs(c).*un.' < 1e-11*norm(fld(Ss))) = 0;
        assert(norm(U*c - tgt) <= 1e-8*norm(fld(Ss)) + 1e-14, 'lambRandomGram: decomposition failed');
        se = zeros(nu, 1);  se(Gr.iOut(~Gr.prop)) = c(~Gr.prop);
        vec = si;  vec(Gr.iOut) = c;
        Gr.states(end+1) = struct('med', mm, 'type', a, 'vec', vec, 'si', si, 'se', se, 'c', c, ...
                                  'N', Ss.I*sqrt(1 - (k*ca/w)^2), 'I', Ss.I);
    end
end
if opts.IW
    SI = lamb.rte.scatteringState(I, 'IW1', [], w, z);
    c = solve(fld(SI));
    assert(norm(U*c - fld(SI)) <= 1e-8*norm(fld(SI)) + 1e-14, 'lambRandomGram: IW decomposition failed');
    Gr.iwVec = zeros(nu, 1);  Gr.iwVec(Gr.iOut) = c;  Gr.iwI = SI.I;
end
% --- outgoing tables at z (as lambEmitLocal / lambStationarity) -----------------------------
kb1 = w/I.lower.beta;  zz = z;  T = M.T;  phys = @(n) n;
if z >= 0
    sd = T.solid;  dh = kb1;
elseif I.upper.isFluid
    sd = T.fluid;  dh = -kb1;
else
    T = M.Tflip;  sd = T.solid;  dh = w/T.iface.lower.beta;  zz = -z;
    s = s.*[1 1 -1 1 1 1 1 -1 -1].';  K = -K;
    phys = @(n) regexprep(regexprep(regexprep(n, '1$', 'X'), '2$', '1'), 'X$', '2');
end
r = w/T.omega0;  sc = [r r r r^2 r^2 r^2 r^2 r^2 r^2];  h = dh*zz;
Kh = K;  if zz < 0, Kh = -K; end                                   % half-line frame y = |z|
nph = opts.Nphi;  phi = 2*pi*(0:nph - 1)/nph;  dphi = 2*pi/nph;
Rm = zeros(9, 9, nph);
for i = 1:nph, Rm(:, :, i) = lamb.rte.rotationMatrix(-phi(i)); end
Bs = permute(reshape(pagemtimes(Rm, reshape(s, 9, [])), 9, nu, 3, nph), [1 4 3 2]);   % 9 x phi x j x u
Gr.G = zeros(nu);
to = {};  if numel(cells) == 2, to = cells(3 - m).types; end
Gr.ch = struct('Q', {}, 'npc', {}, 'prop', {}, 'C', {}, 'name', {}, 'same', {}, 'cell', {});
for C = sd.body
    N = numel(C.k);  nt = size(C.Tr, 2);
    A = C.away.*exp(C.nuAway*h).*sc;  e = C.nuAway;
    for t = 1:nt
        A = cat(3, A, C.Tr(:, t).*reshape(C.tow(:, t, :), N, 9).*exp(-C.nuTow(:, t)*h).*sc);
        e = [e, -C.nuTow(:, t)]; %#ok<AGROW>
    end
    Q = pairForms(A, e, Bs, Kh, dh, zz, C.k*r, k, phi, het).*reshape(C.c, 1, 1, []);
    Q = Q*dphi;
    nm = phys(C.name);  same = any(strcmp(Cm.types, nm));
    if same, cl = find(strcmp(Cm.types, nm)); else, cl = find(strcmp(to, nm)); end
    pr = abs(real(C.nuTow)) <= 1e-9*max(abs(C.nuTow), 1);
    Gr.ch(end+1) = struct('Q', Q, 'npc', 1 + nt, 'prop', pr, 'C', C, 'name', nm, 'same', same, 'cell', cl);
    Gr.G = Gr.G + blockSum(sum(Q, 3), nu);
end
Gr.iwQ = cell(1, numel(sd.iw));
for j = 1:numel(sd.iw)
    P = sd.iw(j);  nt = numel(P.nuTow);  A = zeros(1, 9, nt);
    for t = 1:nt, A(1, :, t) = P.fac*P.q(t)*reshape(P.tow(1, t, :), 1, 9).*exp(-P.nuTow(t)*h).*sc; end
    Gr.iwQ{j} = pairForms(A, -P.nuTow(:).', Bs, Kh, dh, zz, P.k*r, k, phi, het)*dphi;
    Gr.G = Gr.G + blockSum(Gr.iwQ{j}, nu);
end
Gr.G = (Gr.G + Gr.G')/2;
Gr.z = z;  Gr.k = k;  Gr.m = m;  Gr.nu = nu;
end

% ------------------------------------------------------------------------------------------
function Q = pairForms(A, e, Bs, Kh, dh, zz, kap, k, phi, het)
% Q(alpha, beta, n) = sum_{phi, j} T_alpha conj(T_beta) (1/2)[g_alpha + conj g_beta]
persistent hasMex
if isempty(hasMex), hasMex = exist('lambPairFormsMex', 'file') == 3; end   % exist is ~ms: once
if hasMex && ~any(strcmpi(het.law, {'gaussian', 'exp', 'vonkarman'}))     % numerical law: rational g
    [~, R] = lambDepthCorr(het, 0, 0, 0);                          % cached fit (poles, residues)
    Q = lambPairFormsMex(complex(A), complex(e), complex(Bs), complex(Kh(:)), dh, zz, kap(:), k, ...
                         phi(:), complex(R.s(:)), complex(R.r(:)));
    return
end
N = size(A, 1);  npc = size(A, 3);  nph = numel(phi);  nu = size(Bs, 4);  na = nu*npc;
T = zeros(na, nph*3, N);  G = zeros(na, nph, N);
dk = sqrt(max(k^2 + kap(:).^2 - 2*k*kap(:).*cos(phi), 0));         % N x nph
mir = min(0:nph - 1, nph - (0:nph - 1)) + 1;  dk = dk(:, 1:max(mir));   % dk(phi) = dk(2 pi - phi): g on [0, pi]
y = abs(zz);  Kc = -1i*e*dh;  if zz < 0, Kc = -Kc; end            % pieces, half-line frame
for c = 1:npc
    Tc = A(:, :, c)*reshape(Bs, 9, []);                            % N x (nph*3*nu)
    Tc = reshape(Tc, N, nph*3, nu);
    for u = 1:nu
        al = u + nu*(c - 1);
        T(al, :, :) = reshape(Tc(:, :, u).', 1, nph*3, N);
        g = lambDepthCorr(het, Kh(u) + Kc(:, c), y, dk);            % same y, dk for all (c, u): cached
        G(al, :, :) = reshape(g(:, mir).', 1, nph, N);
    end
end
GT = T.*repmat(G, 1, 3, 1);                                        % same g for the 3 columns
Q = (pagemtimes(GT, 'none', T, 'ctranspose') + pagemtimes(T, 'none', GT, 'ctranspose'))/2;
end

function S = blockSum(Q, nu)
% sum over the piece pairs: nu-by-nu
na = size(Q, 1);  npc = na/nu;
S = reshape(sum(sum(reshape(Q, nu, npc, nu, npc), 2), 4), nu, nu);
end

function L = covColumns(med, het)
% C = L L', C_cd = cv_c cv_d r_cd m_c m_d over (lambda, mu, rho)
m = [med.lambda; med.mu; med.rho];  sd = het.cv(:).*m;  r = het.corr;
C = (sd*sd.').*[1 r(1) r(2); r(1) 1 r(3); r(2) r(3) 1];
[V, E] = eig((C + C.')/2);  L = V*sqrt(max(E, 0));
end
