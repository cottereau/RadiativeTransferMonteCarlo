function varargout = lambBulkScatter(varargin)
%LAMBBULKSCATTER  Whole-space Born scattering by point heterogeneities (outside the slab).
%   B = lambBulkScatter('totals', med, het, omega)      per-scatterer cross sections
%   P = lambBulkScatter('scatter', M, P, idx)            scatter the body particles idx
%   D = lambBulkScatter('dsigma', med, het, omega, a, b, ct)   random medium: dsigma/dOmega
%
%   Random medium: het with .Phi (handle of q*a, MaterialClass normalisation
%   Rhat(q) = var a^3 Phi(q a)), .a (correlation length), .cv = [dlam dmu drho]/[lam mu rho]
%   (coefficients of variation), .corr = [r_lm r_lr r_mr]. The point-scatterer n V^2 dm_c dm_d
%   becomes (2 pi)^3 C_cd a^3 Phi(|k_a d - k_b g| a), C_cd = cv_c cv_d r_cd m_c m_d (per unit
%   volume, so rates use n = 1). 'dsigma' gives the phi- and polarisation-averaged differential
%   cross section per unit volume at cos(theta) = ct (incident a = 'P'|'S' along z, out b);
%   'totals' integrates it (adaptive in theta: Phi is not polynomial).
%
%   Equivalent source of a point heterogeneity (V, drho, dlam, dmu) in the unit plane wave
%   u = e e^{i k d.x}:  F = omega^2 drho V e,  M = -V (dlam tr(eps) I + 2 dmu eps),
%   eps = sym(i k d e'). Far field (the same formulas as testLambPhase3 G3.1, where they equal
%   the Lamb whole-space tables to 1e-8):
%       P : A = (g.F - i k_a g'M g)/(4 pi rho alpha^2),   dW/dOmega = rho alpha omega^2 |A|^2/2
%       S : v = (I - g g')(F - i k_b M g)/(4 pi rho beta^2), dW/dOmega = rho beta omega^2 |v|^2/2
%   Totals use an exact product quadrature on the sphere (the integrand is a polynomial of
%   degree <= 4 in g): 8 Gauss nodes in cos(theta) x 16 in phi.
%   'scatter' draws the outgoing type with the precomputed probabilities, the direction by
%   rejection from |A|^2 (bound (|F| + k ||M||_F)^2), and an outgoing S polarisation along
%   e_SV or t (w.r.t. the interface normal) with probabilities |v.e_SV|^2, |v.t|^2, so the
%   observed SV/SH split is exact.
switch varargin{1}
    case 'totals',  varargout{1} = totals(varargin{2:end});
    case 'scatter', varargout{1} = scatter(varargin{2:end});
    case 'dsigma',  varargout{1} = dsigma(varargin{2:end});
end
end

% ------------------------------------------------------------------------------------------
function B = totals(med, het, w)
if isPSD(het)
    types = {'P'};  if ~med.isFluid, types{2} = 'S'; end
    B.sigma = zeros(2);
    for a = 1:numel(types)
        for b = 1:numel(types)
            B.sigma(a, b) = 2*pi*integral(@(ct) dsigma(med, het, w, types{a}, types{b}, ct), -1, 1, ...
                                          'RelTol', 1e-10, 'AbsTol', 0);
        end
    end
    B.pP = B.sigma(:, 1)./sum(B.sigma, 2);
    B.rate = [med.alpha; max(med.beta, eps)].*sum(B.sigma, 2);    % per unit volume: n = 1
    if med.isFluid, B.rate(2) = 0; B.pP(2) = 1; end
    return
end
[x, wx] = gaussLegendre(8);  ph = 2*pi*(0:15)/16;
[CT, PH] = ndgrid(x, ph);  WT = repmat(wx, 1, 16)*(2*pi/16);
g = [sqrt(1 - CT(:).^2).*cos(PH(:)), sqrt(1 - CT(:).^2).*sin(PH(:)), CT(:)];
types = {'P'};  if ~med.isFluid, types{2} = 'S'; end
B.sigma = zeros(2);                                   % rows: incident P, S; cols: out P, S
for a = 1:numel(types)
    d = [0 0 1];  e = d;  if types{a} == 'S', e = [1 0 0]; end
    [F, Mm] = source(med, het, w, types{a}, d, e);
    I = 0.5*med.rho*speed(med, types{a})*w^2;
    for b = 1:numel(types)
        dW = arrayfun(@(i) dWdO(med, w, F, Mm, g(i, :), types{b}), (1:size(g, 1)).');
        B.sigma(a, b) = sum(WT(:).*dW)/I;
    end
end
B.pP = B.sigma(:, 1)./sum(B.sigma, 2);                 % probability of a P outgoing wave
B.rate = het.n*[med.alpha; max(med.beta, eps)].*sum(B.sigma, 2);   % events per unit time
if med.isFluid, B.rate(2) = 0; B.pP(2) = 1; end
end

function P = scatter(M, P, idx)
w = M.omega;
side0 = P.side(idx);  mode0 = P.mode(idx);        % incident type fixed before any conversion
for side = 1:2
    for ty = 'PS'
        sel = idx(side0 == side & mode0 == 1 + (ty == 'S'));
        if isempty(sel), continue, end
        med = M.iface.lower;  if side == 2, med = M.iface.upper; end
        het = M.het{side};  Bk = M.bulk{side};
        n = numel(sel);
        d = P.dir(sel, :);  e = d;  if ty == 'S', e = P.perp(sel, :); end
        if isPSD(het)
            P = scatterPSD(P, sel, med, het, Bk, w, ty, d, e);  continue
        end
        [F, Mm] = source(med, het, w, ty, d, e);                   % F n-by-3, Mm 3-by-3-by-n
        outP = rand(n, 1) < Bk.pP(1 + (ty == 'S'));
        if med.isFluid, outP(:) = true; end
        g = zeros(n, 3);  todo = true(n, 1);
        kout = zeros(n, 1);  kout(outP) = w/med.alpha;  kout(~outP) = w/max(med.beta, eps);
        bound = (vecnorm(F, 2, 2) + kout.*reshape(sqrt(sum(abs(Mm).^2, [1 2])), [], 1)).^2;
        while any(todo)
            t = find(todo);  m = numel(t);
            gg = randn(m, 3);  gg = gg./vecnorm(gg, 2, 2);
            a2 = zeros(m, 1);
            for i = 1:m
                a2(i) = amp2(F(t(i), :), Mm(:, :, t(i)), gg(i, :), kout(t(i)), outP(t(i)));
            end
            ok = rand(m, 1).*bound(t) < a2;
            g(t(ok), :) = gg(ok, :);  todo(t(ok)) = false;
        end
        l = g(:, 1:2);  sl = vecnorm(l, 2, 2);
        l = l./sl;  l(sl == 0, :) = repmat([1 0], nnz(sl == 0), 1);
        tt = [-l(:, 2), l(:, 1), zeros(n, 1)];
        eSV = cross(g, tt, 2);
        pp = cross(g, randn(n, 3), 2);  pp = pp./vecnorm(pp, 2, 2);   % P: any transverse
        for i = find(~outP).'
            v = (eye(3) - g(i, :).'*g(i, :))*(F(i, :).' - 1i*kout(i)*Mm(:, :, i)*g(i, :).');
            if rand*norm(v)^2 < abs(eSV(i, :)*v)^2, pp(i, :) = eSV(i, :); else, pp(i, :) = tt(i, :); end
        end
        P.dir(sel, :) = g;  P.perp(sel, :) = pp;  P.mode(sel) = 1 + ~outP;
    end
end
end

function P = scatterPSD(P, sel, med, het, Bk, w, ty, d, e)
% random medium: outgoing type by Bk.pP, direction by rejection from sum_j |A_j(g)|^2 Shat(q)
% (bound: sum_j bound_j x max Shat on [|k_in - k_out|, k_in + k_out]: Phi need not decrease,
% e.g. the monodisperse law peaks at q Lc ~ 1.25), S polarisation by sum_j |e.v_j|^2
L = covColumns(med, het);  n = numel(sel);  S3 = @(q) (2*pi)^3*het.a^3*het.Phi(het.a*q);
F = cell(1, 3);  Mm = F;
for j = 1:3
    [F{j}, Mm{j}] = source(med, struct('V', 1, 'dlam', L(1, j), 'dmu', L(2, j), 'drho', L(3, j)), w, ty, d, e);
end
kin = w/speed(med, ty);
outP = rand(n, 1) < Bk.pP(1 + (ty == 'S'));  if med.isFluid, outP(:) = true; end
kout = zeros(n, 1);  kout(outP) = w/med.alpha;  kout(~outP) = w/max(med.beta, eps);
bound = zeros(n, 1);
for j = 1:3
    bound = bound + (vecnorm(F{j}, 2, 2) + kout.*reshape(sqrt(sum(abs(Mm{j}).^2, [1 2])), [], 1)).^2;
end
for ko = unique(kout).'                                           % max of Shat over the reachable q
    qq = linspace(abs(kin - ko), kin + ko, 513);
    bound(kout == ko) = bound(kout == ko)*1.05*max(S3(qq));      % ponytail: 513-point grid + 5 % margin
end
g = zeros(n, 3);  todo = true(n, 1);
while any(todo)
    t = find(todo);  m = numel(t);
    gg = randn(m, 3);  gg = gg./vecnorm(gg, 2, 2);
    a2 = zeros(m, 1);
    for i = 1:m
        for j = 1:3, a2(i) = a2(i) + amp2(F{j}(t(i), :), Mm{j}(:, :, t(i)), gg(i, :), kout(t(i)), outP(t(i))); end
        a2(i) = a2(i)*S3(norm(kin*d(t(i), :) - kout(t(i))*gg(i, :)));
    end
    ok = rand(m, 1).*bound(t) < a2;
    g(t(ok), :) = gg(ok, :);  todo(t(ok)) = false;
end
l = g(:, 1:2);  sl = vecnorm(l, 2, 2);
l = l./sl;  l(sl == 0, :) = repmat([1 0], nnz(sl == 0), 1);
tt = [-l(:, 2), l(:, 1), zeros(n, 1)];  eSV = cross(g, tt, 2);
pp = cross(g, randn(n, 3), 2);  pp = pp./vecnorm(pp, 2, 2);
for i = find(~outP).'
    Pg = eye(3) - g(i, :).'*g(i, :);  aSV = 0;  aA = 0;
    for j = 1:3
        v = Pg*(F{j}(i, :).' - 1i*kout(i)*Mm{j}(:, :, i)*g(i, :).');
        aSV = aSV + abs(eSV(i, :)*v)^2;  aA = aA + norm(v)^2;
    end
    if rand*aA < aSV, pp(i, :) = eSV(i, :); else, pp(i, :) = tt(i, :); end
end
P.dir(sel, :) = g;  P.perp(sel, :) = pp;  P.mode(sel) = 1 + ~outP;
end

function L = covColumns(med, het)
% C = L L', C_cd = cv_c cv_d r_cd m_c m_d over (lambda, mu, rho)
m = [med.lambda; med.mu; med.rho];  sd = het.cv(:).*m;
r = het.corr;  C = (sd*sd.').*[1 r(1) r(2); r(1) 1 r(3); r(2) r(3) 1];
[V, E] = eig((C + C.')/2);  L = V*sqrt(max(E, 0));
end

function D = dsigma(med, het, w, a, b, ct)
% phi- and polarisation-averaged dsigma/dOmega per unit volume (random medium), incident along z
L = covColumns(med, het);                                          % C = L L'
d = [0 0 1];  e = d;  if a == 'S', e = [1 0 0]; end               % phi average = pol. average
ph = 2*pi*(0:15)/16;                                               % exact: degree <= 4 in phi
I = 0.5*med.rho*speed(med, a)*w^2;  ka = w/speed(med, a);  kb = w/speed(med, b);
F = cell(1, 3);  Mm = F;
for j = 1:3
    [F{j}, Mm{j}] = source(med, struct('V', 1, 'dlam', L(1, j), 'dmu', L(2, j), 'drho', L(3, j)), w, a, d, e);
end
D = zeros(size(ct));
for i = 1:numel(ct)
    st = sqrt(max(1 - ct(i)^2, 0));  acc = 0;
    for p = ph
        g = [st*cos(p), st*sin(p), ct(i)];
        for j = 1:3, acc = acc + dWdO(med, w, F{j}, Mm{j}, g, b); end
    end
    q = sqrt(max(ka^2 + kb^2 - 2*ka*kb*ct(i), 0));
    D(i) = acc/numel(ph)/I*(2*pi)^3*het.a^3*het.Phi(q*het.a);
end
end

function t = isPSD(het)
t = isfield(het, 'Phi') && ~isempty(het.Phi);
end

function [F, Mm] = source(med, het, w, ty, d, e)
% equivalent source of the unit plane wave (rows of d, e), F n-by-3, Mm 3-by-3-by-n
k = w/speed(med, ty);  n = size(d, 1);
F = w^2*het.drho*het.V*e;
G = 1i*k*reshape(e.', 3, 1, n).*reshape(d.', 1, 3, n);     % du_i/dx_j = i k e_i d_j
epsil = (G + pagetranspose(G))/2;
tr = sum(d.*e, 2)*1i*k;
Mm = -het.V*(het.dlam*reshape(tr, 1, 1, n).*eye(3) + 2*het.dmu*epsil);
end

function a2 = amp2(F, Mm, g, k, isP)
if isP
    a2 = abs(g*F.' - 1i*k*(g*Mm*g.'))^2;
else
    v = (eye(3) - g.'*g)*(F.' - 1i*k*Mm*g.');
    a2 = real(v'*v);
end
end

function dW = dWdO(med, w, F, Mm, g, ty)
c = speed(med, ty);
dW = 0.5*med.rho*c*w^2*amp2(F, Mm, g, w/c, ty == 'P')/(4*pi*med.rho*c^2)^2;
end

function c = speed(med, ty)
if ty == 'P', c = med.alpha; else, c = med.beta; end
end

function [x, w] = gaussLegendre(n)
i = (1:n-1)';  b = i./sqrt(4*i.^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, k] = sort(diag(D));  w = 2*V(1, k)'.^2;
end
