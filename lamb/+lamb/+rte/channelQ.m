function [Q, W] = channelQ(T, zeta, omega, s)
%CHANNELQ  Hermitian power forms of every channel: W_c = s' * Q_c * s.
%   Q = lamb.rte.channelQ(T, zeta, omega)
%   [Q, W] = lamb.rte.channelQ(T, zeta, omega, s)
%
%   Q.body.(name)  9-by-9 for each body channel (P1, SV1, SH1, P2, SV2, SH2, Pf)
%   Q.iw{j}        9-by-9 for each real interface-wave pole
%   Q.leaky{j}     9-by-9: part of the body forms inside the Lorentzian peak of near-real
%                  leaky pole j (already included in Q.body; quasi-trapped energy)
%   W              (with s) the powers: .body.(name) .iw(j) .leaky(j) .bodyTotal .ifTotal
%                  .total (= bodyTotal + ifTotal) and .pct (same layout, % of total)
%
%   The azimuth integral uses 9 uniform nodes, which is exact: the amplitudes are
%   trigonometric polynomials of degree <= 2 in phi, so |.|^2 has degree <= 4.
Nphi = 9;  phi = 2*pi*(0:Nphi-1)/Nphi;
Rm = arrayfun(@(p) lamb.rte.rotationMatrix(-p), phi, 'UniformOutput', false);
avg = @(A) azimuthAverage(A, Rm, 2*pi/Nphi);

sd = T.solid;
if zeta < 0, sd = T.fluid; end
Q.leaky = repmat({zeros(9)}, 1, numel(T.leaky));
for C = sd.body
    [d, c] = lamb.rte.channelAmplitudes(T, C.name, zeta, omega);
    Q.body.(C.name) = avg(d'*(c.*d));
    for j = 1:numel(T.leaky)
        m = C.leakMask(:, j);
        Q.leaky{j} = Q.leaky{j} + avg(d(m, :)'*(c(m).*d(m, :)));
    end
end
Q.iw = cell(1, numel(sd.iw));
for j = 1:numel(sd.iw)
    hj = lamb.rte.channelAmplitudes(T, sprintf('IW%d', j), zeta, omega);
    Q.iw{j} = avg(hj'*hj);
end

if nargin == 4
    f = @(A) real(s'*A*s);
    W.body  = structfun(f, Q.body, 'UniformOutput', false);
    W.iw    = cellfun(f, Q.iw);
    W.leaky = cellfun(f, Q.leaky);
    W.bodyTotal = sum(cell2mat(struct2cell(W.body)));
    W.ifTotal   = sum(W.iw);
    W.total     = W.bodyTotal + W.ifTotal;
    p = @(v) 100*v/W.total;
    W.pct.body  = structfun(p, W.body);
    W.pct.body  = cell2struct(num2cell(W.pct.body), fieldnames(W.body));
    W.pct.iw    = p(W.iw);
    W.pct.leaky = p(W.leaky);
    W.pct.trapped = p(W.ifTotal + sum(W.leaky));
end
end

function Q = azimuthAverage(A, Rm, dphi)
% int dphi R(-phi).' * A * R(-phi), exact for the degree-4 trigonometric content
Q = zeros(9);
for i = 1:numel(Rm)
    Q = Q + Rm{i}.'*A*Rm{i};
end
Q = dphi*Q;
Q = (Q + Q')/2;                                       % Hermitian to rounding
end
