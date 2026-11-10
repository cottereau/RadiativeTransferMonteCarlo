function [g, R] = lambDepthCorr(het, K, y, dk)
%LAMBDEPTHCORR  One-sided depth correlation of a random medium cut at the interface.
%   g = lambDepthCorr(het, K, y, dk)
%
%   g(K, y; dk) = int_{-y}^{inf} e^{i K zeta} C(dk, zeta) dzeta,
%   C(dk, zeta) = (1/2pi) int dq Shat(sqrt(dk^2 + q^2)) e^{i q zeta},  Shat = (2 pi)^3 a^3 Phi(a q)
%   (MaterialClass normalisation, the covariance is applied elsewhere): the weight that a term
%   ~ e^{i K z'} of the field at depth z' in (0, inf) gets at the depth y >= 0 through the
%   correlation (lambRandomKernels, 'Qz','nonlocal'). Point limit (a -> 0): Shat(0), constant;
%   far from the interface (y >> a): Shat at the complex q = (dk, K).
%   Closed forms for het.law = 'gaussian' and 'exp' ('VonKarman' with het.nu = 0.5 is 'exp');
%   any other het.law (e.g. 'monodispersesphere', 'Imported'): numerical, for any het.Phi.
%   Im K >= 0 (terms that do not grow away from the interface). K, y, dk broadcast.
%
%   Numerical law: one rational fit F(s) = Shat(sqrt(s)) ~ sum_l r_l/(s - s_l) on s = q^2 >= 0
%   (AAA, poles off [0, inf), residues by least squares, |error| < 1e-5 max F; see rationalFit for
%   the two sample sets). Then Shat(sqrt(dk^2 + q^2)) = sum r_l/(q^2 + P_l^2),
%   P_l^2 = dk^2 - s_l, so C = sum r_l e^{-P_l |zeta|}/(2 P_l) and
%   g = sum r_l/(2 P_l) [1/(P_l - iK) + (1 - e^{-(P_l + iK) y})/(P_l + iK)]   (closed form, any dk).
%   The fit is cached (persistent) per (a, Phi).
%
%   References:
%     - Nakatsukasa, Y., Sete, O. & Trefethen, L.N. (2018), The AAA algorithm for rational
%       approximation, SIAM J. Sci. Comput. 40, A1494-A1522. doi:10.1137/16M1106122
%     - Weideman, J.A.C. (1994), Computation of the complex error function, SIAM J. Numer.
%       Anal. 31, 1497-1518. doi:10.1137/0731077
switch lower(het.law)
    case 'gaussian'
        ag = (3*sqrt(pi)/4)^(2/3);  b = het.a^2/(4*ag);             % Shat = F0 e^{-b (dk^2 + q^2)}
        F0 = (2*pi)^3*het.a^3/(6*pi^2);
        xi = sqrt(b)*K - 1i*y/(2*sqrt(b));                          % gamma = int_{-y}^inf e^{iK zeta} N(0, 2b)
        gam = exp(-b*K.^2) - 0.5*exp(-y.^2/(4*b) - 1i*K.*y).*faddeeva(-xi);
        g = F0*exp(-b*dk.^2).*gam;
    case {'exp', 'vonkarman'}
        if strcmpi(het.law, 'vonkarman') && abs(het.nu - 0.5) > 1e-12
            error('lambDepthCorr:nu', 'von Karman slab media: only nu = 0.5 (= exponential)');
        end
        al = 6^(1/3);  F0 = (2*pi)^3*het.a^3/(6*pi^2)*(al/het.a)^4;   % Shat = F0/(s^2 + q^2)^2
        s = sqrt(al^2/het.a^2 + dk.^2);                             % C = F0 (1 + s|z|) e^{-s|z|}/(4 s^3)
        p = s - 1i*K;  r = s + 1i*K;
        Ip = 1./p + s./p.^2;                                        % zeta > 0
        rz = r.*y;  Im = (1 - exp(-rz))./r + s.*(1 - exp(-rz).*(1 + rz))./r.^2;   % 0 > zeta > -y
        sm = abs(rz) < 1e-3;                                        % series where r y is tiny
        if any(sm(:))
            [rr, ss, yy] = deal(r + 0*rz, s + 0*rz, y + 0*rz);
            rr = rr(sm);  ss = ss(sm);  yy = yy(sm);
            Im(sm) = yy + (ss - rr).*yy.^2/2 + (rr.^2/2 - ss.*rr).*yy.^3/3 + (-rr.^3/6 + ss.*rr.^2/2).*yy.^4/4;
        end
        g = F0*(Ip + Im)./(4*s.^3);
    otherwise
        R = rationalFit(het);  g = rationalG(R, K, y, dk);          % R: the fit (poles, residues)
end
end

% ------------------------------------------------------------------------------------------
function R = rationalFit(het)
% poles s_l and residues r_l of F(s) = Shat(sqrt(s)), cached for the last (a, Phi) seen.
% Two fits: 'full' (dense to q Lc 60, true samples to 1e4) and 'cut' (dense up to where the envelope
% of |F| < 1e-6 max, then F = 0 anchors: an oscillating tail sampled sparsely costs poles). g needs
% F continued to complex q, which the cut fit can spoil (e.g. Gaussian): the cut fit is used only
% if its g matches the full fit's within 1e-6 over K in [-4, 4] + i [0, 3], y in [0, 10], dk in [0, 4].
persistent cache                                                    % called per event: keep the hit cheap
if ~isempty(cache) && cache.a == het.a && isequal(cache.Phi, het.Phi), R = cache.R; return, end
a = het.a;  F = @(s) (2*pi)^3*a^3*het.Phi(a*sqrt(s));
x = linspace(0, 200, 20001);  f = F((x/a).^2);                      % q Lc samples, dense
if any(~isfinite(f)), error('lambDepthCorr:Phi', 'het.Phi is not finite on [0, %g]', x(end)); end
xt = [x, logspace(log10(x(end)), 4, 400)];  ft = F((xt/a).^2);     % for the error checks
R = fitOn(xt(xt <= 60 | xt > 200), ft(xt <= 60 | xt > 200), a, xt, ft);
n = find(cummax(abs(f), 'reverse') >= 1e-6*max(abs(f)), 1, 'last') + 1;
if n < numel(x)
    xz = logspace(log10(x(n)*1.01), 4, 200);
    Rc = fitOn([x(1:n), xz], [f(1:n), zeros(size(xz))], a, xt, ft);
    [Kr, Ki, yy, dk] = ndgrid(-4:0.25:4, [0 0.5 1 3], [0 0.3 1 3 10]*a, 0:0.5:4);
    gf = rationalG(R, Kr + 1i*Ki, yy, dk);  gc = rationalG(Rc, Kr + 1i*Ki, yy, dk);
    if max(abs(gc(:) - gf(:))) < 1e-6*max(abs(gf(:))), R = Rc; end
end
cache = struct('a', a, 'Phi', het.Phi, 'R', R);
end

function R = fitOn(x, f, a, xt, ft)
% AAA in t = 1/(1 + s a^2) through F(inf) = 0, residues by least squares; checked on (xt, ft)
s = (x/a).^2;  t = [0, 1./(1 + s*a^2)];  fa = [0, f];
[w, zt] = aaa(fa(:), t(:), 1e-8*max(abs(f)), 400);
pt = barPoles(w, zt);  sl = (1./pt - 1)/a^2;                         % poles in t, then in s
sl = sl(isfinite(sl) & ~(abs(imag(sl)) < 1e-9*max(1, abs(sl)) & real(sl) >= 0));
r = (1./(s(:) - sl.'))\f(:);
st = (xt/a).^2;  err = max(abs(real(sum(r.'./(st(:) - sl.'), 2)) - ft(:)))/max(abs(ft));
if err > 1e-5
    error('lambDepthCorr:fit', 'rational fit of Phi: max error %.1e of max Shat (> 1e-5)', err);
end
R = struct('s', sl, 'r', r, 'err', err);
end

function g = rationalG(R, K, y, dk)
% P_l = sqrt(dk^2 - s_l), r_l/(2 P_l) and e^{-P_l y} depend on (dk, y) only: kept for the next
% call with the same (dk, y) (lambRandomGram calls with one (dk, y) and many K).
% 1/(P - iK) + (1 - e^{-(P + iK) y})/(P + iK) = [P + iK + (1 - e^{-Py} e^{-iKy}) (P - iK)]/(P^2 + K^2)
persistent pre hasMex
if isempty(hasMex), hasMex = exist('lambRationalGmex', 'file') == 3; end   % exist is ~ms: once
sz = size(K + y + dk);
if isempty(pre) || ~isequal(pre.s, R.s) || ~isequal(pre.dk, dk)
    P2 = dk(:).^2 - R.s.';  P = sqrt(P2);                           % Re P > 0: s_l off [0, inf)
    pre = struct('s', R.s, 'dk', dk, 'P', P, 'P2', P2, 'A', R.r.'./(2*P), 'y', NaN, 'E', [], 'Et', []);
    pre.Pt = complex(P.');  pre.P2t = complex(P2.');  pre.At = complex(pre.A.');   % poles contiguous (MEX)
end
if isscalar(y) && ~isequal(pre.y, y), pre.y = y;  pre.E = exp(-pre.P*y);  pre.Et = complex(pre.E.'); end
K = reshape(K + zeros(sz), [], 1);  y = reshape(y + zeros(sz), [], 1);
whole = isequal(sz, size(dk)) || numel(dk) == 1;                    % dk entry q is query q (or all)
if ~whole, id = reshape(reshape(1:numel(dk), size(dk)) + zeros(sz), [], 1); end
if hasMex                                                            % compiled inner loop
    if ~whole, jd = id; elseif numel(dk) > 1, jd = (1:numel(K)).'; else, jd = ones(numel(K), 1); end
    Et = [];  if ~isempty(pre.Et) && all(y == pre.y), Et = pre.Et; end
    g = reshape(lambRationalGmex(pre.Pt, pre.P2t, pre.At, Et, jd, complex(K), y), sz);
    return
end
nc = max(1, floor(2^18/numel(R.s)));                                % queries per chunk
g = complex(zeros(numel(K), 1));
for c0 = 1:nc:numel(K)
    ii = c0:min(c0 + nc - 1, numel(K));
    if whole && numel(dk) > 1, jj = ii; elseif whole, jj = ones(size(ii)); else, jj = id(ii); end
    if numel(ii) == numel(K) && whole && numel(dk) > 1               % one chunk: no copies
        P = pre.P;  P2 = pre.P2;  A = pre.A;  E = pre.E;
    else
        P = pre.P(jj, :);  P2 = pre.P2(jj, :);  A = pre.A(jj, :);  E = [];
        if ~isempty(pre.E), E = pre.E(jj, :); end
    end
    yi = y(ii);  iK = 1i*K(ii);
    if isempty(E) || ~all(yi == pre.y), E = exp(-P.*yi); end
    t = P + iK;
    v = A.*(t + (1 - E.*exp(-iK.*yi)).*(P - iK))./(P2 + K(ii).^2);
    sm = abs(t).*yi < 1e-4;                                         % small (P + iK) y: series
    if any(sm(:))
        x = t.*yi;  yy = yi + 0*x;  Pm = P - iK + 0*x;  Aa = A + 0*x;
        v(sm) = Aa(sm).*(1./Pm(sm) + yy(sm).*(1 - x(sm)/2 + x(sm).^2/6));
    end
    g(ii) = sum(v, 2);
end
g = reshape(g, sz);
end

function [w, z, fz] = aaa(F, Z, tol, mmax)
% AAA rational approximation (Nakatsukasa, Sete & Trefethen 2018); support point 1 forced
M = numel(Z);  J = (1:M).';  z = zeros(0, 1);  fz = z;  C = zeros(M, 0);  R = zeros(M, 1);
for m = 1:mmax
    if m == 1, j = 1; else, [~, j] = max(abs(F - R)); end
    z(m, 1) = Z(j);  fz(m, 1) = F(j);  J(J == j) = [];
    C(:, m) = 1./(Z - Z(j));
    A = F.*C - C.*fz.';
    [~, ~, V] = svd(A(J, :), 0);  w = V(:, m);
    N = C*(w.*fz);  D = C*w;  R = F;  R(J) = N(J)./D(J);
    if max(abs(F - R)) <= tol, break, end
end
end

function p = barPoles(w, z)
m = numel(w);  B = eye(m + 1);  B(1, 1) = 0;
E = [0, w.'; ones(m, 1), diag(z)];
p = eig(E, B);  p = p(isfinite(p));
end

function w = faddeeva(z)
% w(z) = e^{-z^2} erfc(-i z) for Im z >= 0 (Weideman 1994, N = 32); Im z < 0 by reflection
w = zeros(size(z));  lo = imag(z) < 0;
w(~lo) = cef(z(~lo));
if any(lo(:)), w(lo) = 2*exp(-z(lo).^2) - cef(-z(lo)); end
end

function w = cef(z)
N = 32;  M = 2*N;  M2 = 2*M;  k = (-M + 1:M - 1).';
L = sqrt(N/sqrt(2));  t = L*tan(k*pi/M/2);
f = [0; exp(-t.^2).*(L^2 + t.^2)];
a = real(fft(fftshift(f)))/M2;  a = flipud(a(2:N + 1));
Z = (L + 1i*z)./(L - 1i*z);
w = 2*polyval(a, Z)./(L - 1i*z).^2 + (1/sqrt(pi))./(L - 1i*z);
end
