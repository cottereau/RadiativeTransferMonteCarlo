function K = lambRandomKernels(T, inc, p, omega, z, het, opts)
%LAMBRANDOMKERNELS  Near-interface scattering of a random medium (PSD), per unit volume.
%   K = lambRandomKernels(T, inc, p, omega, z, het)
%   K = lambRandomKernels(T, inc, p, omega, z, het, 'Tflip', Tf, 'Nphi', 64, 'Qz', 'real')
%
%   Random-medium counterpart of lamb.rte.slabKernels (same T, inc, p, omega, z, Tflip; same
%   outputs K.sigma.(c), K.iw, K.total, K.I), with het the random medium of the depth z as in
%   lambBulkScatter (.Phi .a .cv .corr). Cross sections are per unit volume (point-scatterer
%   n sigma).
%
%   The state at z is a sum of plane-wave components n (incident wave, outgoing waves of the
%   medium of z, evanescent included), each with its equivalent source s_n and vertical
%   wavenumber K_n (s_n ~ e^{i K_n z}); the outgoing amplitude row of a channel node is a sum of
%   pieces c (direct, via each toward wave t), A_c ~ e^{i K_c z}. For a point scatterer the
%   pattern is |sum_{n,c} A_c R s_n|^2; for a stationary random field every pair of terms
%   alpha = (n,c), beta = (n',c') is weighted by the spectrum at the pair's wavevector transfer,
%       (2 pi)^3 a^3 Phi(a |q|),  |q|^2 = |k x - kappa e_phi|^2 + q_z^2,
%       q_z = (K_alpha + conj(K_beta))/2              (complex for evanescent terms; 'full')
%   or q_z = Re(...) ('real': the evanescent decay rates dropped, the interfaceFormFactors
%   approximation, q = k_b n - k_I x with a horizontal guided wavevector; gate G5.3). The pair form
%   assumes the correlation length is short against the distance to the interface (both points
%   of the correlation in the medium of z). As a -> 0 it is the point scatterer with
%   n V^2 dm_c dm_d = (2 pi)^3 a^3 Phi(0) C_cd (gate G5.1).
%
%   'Depth','integrated': rates integrated over the whole medium of z (z only picks the side;
%   every term must decay away from the interface, i.e. an interface wave in or out; channels
%   with a non-decaying term, body -> body, come back NaN). With
%   X_alpha = X_alpha(0) e^{i K_alpha z} on the half line, the double depth integral is
%       (1/2 pi) int dq_z Shat(q) / ((q_z - K_alpha)(q_z - conj K_beta))     ('Qz','exact')
%   (real q_z axis, no full-line assumption); 'full' / 'real' give the local form integrated in
%   closed form, Shat(q_mid) i/(K_alpha - conj K_beta). Gate G5.3.
%
%   'Qz','nonlocal' (per depth; het.law 'gaussian' | 'exp'): the pair weight is
%   (1/2)[g(K_alpha, y) + conj g(K_beta, y)], g the one-sided depth correlation cut at the
%   interface (lambDepthCorr), y = |z|. Its depth integral is the 'exact' rate; it is the point
%   scatterer as a -> 0 and the local 'full' form for y >> a.
arguments
    T struct
    inc char
    p double
    omega (1,1) double
    z double
    het struct
    opts.Tflip = []
    opts.Nphi (1,1) double = 64
    opts.Qz char {mustBeMember(opts.Qz, {'full', 'real', 'exact', 'nonlocal'})} = 'full'
    opts.Depth char {mustBeMember(opts.Depth, {'local', 'integrated'})} = 'local'
    opts.Nq (1,1) double = 256
end
if strcmp(opts.Qz, 'exact') && strcmp(opts.Depth, 'local')
    error('lambRandomKernels:exact', '''Qz'',''exact'' needs ''Depth'',''integrated''.');
end
if strcmp(opts.Qz, 'nonlocal') && strcmp(opts.Depth, 'integrated')
    error('lambRandomKernels:nonlocal', '''Qz'',''nonlocal'' is per depth (its integral is ''exact'').');
end
md = struct('qz', opts.Qz, 'int', strcmp(opts.Depth, 'integrated'), 'Nq', opts.Nq);
I = T.iface;  z = z(:).';  w = omega;
kb1 = w/I.lower.beta;
phi = 2*pi*(0:opts.Nphi - 1)/opts.Nphi;  dphi = 2*pi/opts.Nphi;
Rm = zeros(9, 9, opts.Nphi);
for i = 1:opts.Nphi, Rm(:, :, i) = lamb.rte.rotationMatrix(-phi(i)); end
ss = strcmp(I.kind, 'solid-solid');
S0 = lamb.rte.scatteringState(I, inc, p, w, z(1));
K.z = z;  K.I = S0.I;  K.iw = zeros(0, numel(z));  K.leaky = [];
for iz = 1:numel(z)
    [s, Kn, k] = components(I, inc, p, w, z(iz), het);            % 9 x nc x 3, 1 x nc
    zz = z(iz);  Tz = T;  names = @(n) n;
    if zz >= 0
        sd = T.solid;  dh = kb1;
    elseif ~ss
        sd = T.fluid;  dh = -kb1;
    else                                                           % upper solid: mirrored problem
        if isempty(opts.Tflip), error('lambRandomKernels:flip', 'Scatterers in an upper solid need Tflip.'); end
        Tz = opts.Tflip;  sd = Tz.solid;  dh = w/Tz.iface.lower.beta;  zz = -zz;
        s = s.*[1 1 -1 1 1 1 1 -1 -1].';  Kn = -Kn;
        names = @(n) regexprep(regexprep(regexprep(n, '1$', 'X'), '2$', '1'), 'X$', '2');
    end
    r = w/Tz.omega0;  sc = [r r r r^2 r^2 r^2 r^2 r^2 r^2];  h = dh*zz;
    B = permute(pagemtimes(Rm, reshape(s, 9, [])), [1 3 2]);      % 9 x Nphi x (nc*3)
    for C = sd.body
        N = numel(C.k);  nt = size(C.Tr, 2);
        A = C.away.*exp(C.nuAway*h).*sc;  e = C.nuAway;           % pieces: N x 9 x npc
        for t = 1:nt
            A = cat(3, A, C.Tr(:, t).*reshape(C.tow(:, t, :), N, 9).*exp(-C.nuTow(:, t)*h).*sc);
            e = [e, -C.nuTow(:, t)]; %#ok<AGROW>
        end
        W = pattern(A, e, B, Kn, dh, zz, C.k*r, k, phi, het, md);   % N x Nphi
        K.sigma.(names(C.name))(iz) = sum(C.c.*W, 'all')*dphi/K.I;
    end
    for j = 1:numel(sd.iw)
        P = sd.iw(j);  nt = numel(P.nuTow);  A = zeros(1, 9, nt);
        for t = 1:nt, A(1, :, t) = P.fac*P.q(t)*reshape(P.tow(1, t, :), 1, 9).*exp(-P.nuTow(t)*h).*sc; end
        W = pattern(A, -P.nuTow(:).', B, Kn, dh, zz, P.k*r, k, phi, het, md);
        K.iw(j, iz) = sum(W)*dphi/K.I;
    end
end
K.total = sum(cell2mat(struct2cell(K.sigma)), 1) + sum(K.iw, 1);
end

% ------------------------------------------------------------------------------------------
function W = pattern(A, e, B, Kn, dh, zz, kap, k, phi, het, md)
% sum over pairs of X_alpha conj(X_beta) (2 pi)^3 a^3 Phi(a q); A N x 9 x npc rows with
% exponent coefficients e (N x npc, A ~ e^{e h}); B 9 x Nphi x (nc*3) rotated sources
N = size(A, 1);  npc = size(A, 3);  nph = size(B, 2);  nc = numel(Kn);
X = zeros(N, nph, nc, 3, npc);
for c = 1:npc
    X(:, :, :, :, c) = reshape(A(:, :, c)*reshape(B, 9, []), N, nph, nc, 3);
end
X = reshape(permute(X, [1 2 4 3 5]), N, nph, 3, nc*npc);          % alpha = n + nc (c - 1)
Ka = reshape(Kn(:).' + zeros(N, 1), N, nc, 1) + reshape(-1i*e*dh, N, 1, npc);   % N x nc x npc
Ka = reshape(Ka, N, nc*npc);
if md.int                                                          % X at the interface, y = |z|
    X = X.*reshape(exp(-1i*Ka*zz), N, 1, 1, []);
    if zz < 0, Ka = -Ka; end
    if any(imag(Ka(:)) <= 0), W = nan(N, nph); return, end        % body -> body: infinite
end
dk2 = k^2 + kap(:).^2 - 2*k*kap(:).*cos(phi);                       % N x Nphi
W = zeros(N, nph);  na = nc*npc;  f = (2*pi)^3*het.a^3;
if strcmp(md.qz, 'nonlocal')                                       % one-sided, symmetrised
    Kh = Ka;  if zz < 0, Kh = -Ka; end                             % half-line frame y = |z|
    dk = sqrt(max(dk2, 0));  G = zeros(N, nph, na);
    for al = 1:na, G(:, :, al) = lambDepthCorr(het, Kh(:, al), abs(zz), dk); end
    for al = 1:na
        for be = 1:na
            xx = sum(X(:, :, :, al).*conj(X(:, :, :, be)), 3);
            W = W + real(xx.*(G(:, :, al) + conj(G(:, :, be)))/2);
        end
    end
    return
end
if md.int && strcmp(md.qz, 'exact')
    % 1/((q-Ka)(q-conj Kb)) = [1/(q-Ka) - 1/(q-conj Kb)]/(Ka - conj Kb); single-pole integrals
    % J(K) = (1/2pi) int_R Shat/(q-K) = (1/2pi) int_0^inf Shat 2K/(q^2-K^2) (Shat even in q)
    Kall = reshape([Ka, conj(Ka)], N, 1, 2*na);
    J = integral(@(q) f*het.Phi(het.a*sqrt(dk2 + q^2)).*(2*Kall./(q^2 - Kall.^2)), 0, Inf, ...
                 'ArrayValued', true, 'RelTol', 1e-9, 'AbsTol', 0)/(2*pi);
    for al = 1:na
        for be = 1:na
            xx = sum(X(:, :, :, al).*conj(X(:, :, :, be)), 3);
            W = W + real(xx.*(J(:, :, al) - J(:, :, na + be))./(Ka(:, al) - conj(Ka(:, be))));
        end
    end
    return
end
for al = 1:na
    for be = 1:na
        xx = sum(X(:, :, :, al).*conj(X(:, :, :, be)), 3);           % N x Nphi (sum over cov columns)
        qz_ = (Ka(:, al) + conj(Ka(:, be)))/2;  if strcmp(md.qz, 'real'), qz_ = real(qz_); end
        q2 = dk2 + qz_.^2;
        g = f.*het.Phi(het.a*sqrt(q2));
        if md.int, g = g.*(1i./(Ka(:, al) - conj(Ka(:, be)))); end       % int_0^inf dZ e^{i(Ka - conj Kb)Z}
        W = W + real(xx.*g);
    end
end
end

function [s, Kn, k] = components(I, inc, p, w, z, het)
% plane-wave components of the scattering state at z (medium of z): equivalent sources of the
% three covariance columns (9 x nc x 3) and vertical wavenumbers Kn (s_n ~ e^{i Kn z})
cells = lambMediumCells(I);  m = 1 + (z < 0);  sm = 3 - 2*m;       % +1 lower, -1 upper
med = I.lower;  if m == 2, med = I.upper; end
L = covColumns(med, het);
St = lamb.rte.scatteringState(I, inc, p, w, z);  k = St.k;  pk = k/w;
fld = @(X) [X.u(:); X.eps(:)];
Cm = cells(m);  nq = numel(Cm.types);
U = zeros(12, nq);  waves = cell(1, nq);  Kn = zeros(1, 0);
for q = 1:nq
    waves{q} = lamb.rte.scatteringState(I, Cm.types{q}, pk, w, z, 'Mode', 'incident', 'Away', true);
    U(:, q) = fld(waves{q});
end
tgt = fld(St);  parts = {};
if any(strcmp(Cm.types, inc))                                      % incident wave of this medium
    Si = lamb.rte.scatteringState(I, inc, p, w, z, 'Mode', 'incident');
    tgt = tgt - fld(Si);  parts{end+1} = Si;
    ci = Cm.c(strcmp(Cm.types, inc));  Kn(end+1) = -sm*sqrt(complex(w^2/ci^2 - k^2));
end
un = vecnorm(U);  c = ((U./un)\tgt)./un.';
assert(norm(U*c - tgt) <= 1e-8*norm(fld(St)) + 1e-14, 'lambRandomKernels: decomposition failed');
for q = 1:nq
    if abs(c(q))*un(q) < 1e-11*norm(fld(St)), continue, end
    Sq = waves{q};  Sq.u = c(q)*Sq.u;  Sq.eps = c(q)*Sq.eps;  parts{end+1} = Sq; %#ok<AGROW>
    Kn(end+1) = sm*sqrt(complex(w^2/Cm.c(q)^2 - k^2)); %#ok<AGROW>
end
s = zeros(9, numel(parts), 3);
for n = 1:numel(parts)
    for j = 1:3
        hj = struct('V', 1, 'n', 1, 'dlam', L(1, j), 'dmu', L(2, j), 'drho', L(3, j));
        s(:, n, j) = lamb.rte.equivalentSource(parts{n}, hj, w);
    end
end
end

function L = covColumns(med, het)
% C = L L', C_cd = cv_c cv_d r_cd m_c m_d over (lambda, mu, rho)
m = [med.lambda; med.mu; med.rho];  sd = het.cv(:).*m;  r = het.corr;
C = (sd*sd.').*[1 r(1) r(2); r(1) 1 r(3); r(2) r(3) 1];
[V, E] = eig((C + C.')/2);  L = V*sqrt(max(E, 0));
end
