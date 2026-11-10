function [d, c, k, dOm] = channelAmplitudes(T, name, zeta, omega)
%CHANNELAMPLITUDES  Radiated amplitudes of the 9 basis sources in one channel, azimuth 0.
%   [d, c, k, dOm] = lamb.rte.channelAmplitudes(T, name, zeta, omega)
%
%   T      tables from lamb.rte.radiationTables
%   name   body channel ('P1' 'SV1' 'SH1' 'P2' 'SV2' 'SH2' 'Pf') or 'IW1', 'IW2', ...
%   zeta   k_beta1 * z of the source (>= 0 solid side, < 0 fluid side)
%   omega  frequency (any; the tables scale from T.omega0)
%
%   Body channel: d is N-by-9 (one row per quadrature node k(n)); for a source s the power is
%       W = int dphi sum_n c(n) |d(n,:) * rotationMatrix(-phi) * s|^2,
%   c(n) = w_n k_n E_m(k_n)/(4 pi^2) holding the node weight, so c(n)|d(n,:) R(-phi) s|^2 is
%   the discrete density for sampling (k_n, phi). k is in physical units at omega.
%   dOm(n) is the per-solid-angle weight: dW/dOmega in the direction of node n and azimuth phi
%   is dOm(n)|d(n,:) R(-phi) s|^2, consistent at any omega (c and dOm are in table units and
%   all the frequency scaling is in d).
%   Interface wave: d is 1-by-9 with dW/dphi = |d * R(-phi) * s|^2; c = 1, k = k_p.
sd = T.solid;  h = zeta;
if zeta < 0, sd = T.fluid;  h = -zeta; end
r  = omega/T.omega0;
sc = [r r r r^2 r^2 r^2 r^2 r^2 r^2];                 % force ~ omega, moment ~ omega^2
if strncmp(name, 'IW', 2)
    P = sd.iw(str2double(name(3:end)));
    d = P.fac*reshape(sum(P.q.*P.tow.*exp(-P.nuTow*h), 2), 1, 9).*sc;
    c = 1;  k = P.k*r;  dOm = 1;
    return
end
C = sd.body(strcmp({sd.body.name}, name));
if isempty(C), error('lamb:channelAmplitudes:name', 'No channel %s on this side.', name); end
d = reshape(sum(C.Tr.*C.tow.*exp(-C.nuTow*h), 2), [], 9) + C.away.*exp(C.nuAway*h);
d = d.*sc;
c = C.c;  k = C.k*r;  dOm = C.dOmega;
end
