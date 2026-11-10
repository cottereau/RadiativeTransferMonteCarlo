function [W, Wc] = totalPower(T, zeta, omega, s)
%TOTALPOWER  Radiated power of many point sources at many depths, vectorised (one side).
%   [W, Wc] = lamb.rte.totalPower(T, zeta, omega, s)
%
%   zeta  1-by-Ne source depths k_beta1*z, all on the same side (all >= 0 or all < 0)
%   s     9-by-Ne source vectors (lamb.rte.sourceVector), one per depth
%   W     1-by-Ne total power (all body channels + real interface waves), = sum of the channel
%         powers of lamb.rte.channelQ for each source
%   Wc    struct: .body.(name) 1-by-Ne, .iw nIW-by-Ne, .leaky nL-by-Ne (part of .body.Pf)
%   Same exact 9-node azimuth average as channelQ; the depth dependence is evaluated
%   analytically for all depths at once (pagemtimes over depths).
zeta = zeta(:).';  Ne = numel(zeta);
sd = T.solid;  h = zeta;
if all(zeta < 0), sd = T.fluid;  h = -zeta; end
r  = omega/T.omega0;
sc = [r r r r^2 r^2 r^2 r^2 r^2 r^2].';
s  = s.*sc;                                           % frequency scaling moved onto s
phi = 2*pi*(0:8)/9;
A = zeros(9, 9, Ne);                                  % A(:, m, e) = R(-phi_m) s_e
for m = 1:9
    A(:, m, :) = reshape(lamb.rte.rotationMatrix(-phi(m))*s, 9, 1, Ne);
end
W = zeros(1, Ne);
Wc.leaky = zeros(numel(T.leaky), Ne);
for C = sd.body
    N = numel(C.k);
    d = reshape(C.away, N, 9, 1).*reshape(exp(C.nuAway*h), N, 1, Ne);
    for t = 1:size(C.Tr, 2)
        d = d + (C.Tr(:, t).*reshape(C.tow(:, t, :), N, 9)).*reshape(exp(-C.nuTow(:, t)*h), N, 1, Ne);
    end
    P = abs(pagemtimes(d, A)).^2;                     % N-by-9-by-Ne
    pw = reshape(sum(C.c.*sum(P, 2), 1), 1, Ne)*(2*pi/9);
    Wc.body.(C.name) = pw;
    W = W + pw;
    for j = 1:numel(T.leaky)
        m = C.leakMask(:, j);
        Wc.leaky(j, :) = Wc.leaky(j, :) + reshape(sum(C.c(m).*sum(P(m, :, :), 2), 1), 1, Ne)*(2*pi/9);
    end
end
Wc.iw = zeros(numel(sd.iw), Ne);
for j = 1:numel(sd.iw)
    Pj = sd.iw(j);
    hj = Pj.fac*sum(Pj.q.*Pj.tow.*reshape(exp(-Pj.nuTow.'*h), 1, [], 1, Ne), 2);  % 1-by-1-by-9-by-Ne
    hj = reshape(hj, 1, 9, Ne);
    Wc.iw(j, :) = reshape(sum(abs(pagemtimes(hj, A)).^2, 2), 1, Ne)*(2*pi/9);
    W = W + Wc.iw(j, :);
end
end
