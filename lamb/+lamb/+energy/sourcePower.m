function Ws = sourcePower(iface, src, omega, W0, Wi, poles, leaky, opts)
%SOURCEPOWER  Work done by the source, independent of the flux bookkeeping (eq:Wsrc3).
%   Ws = lamb.energy.sourcePower(iface, src, omega, W0, Wi, poles, leaky, opts)
%
%   Ws.total = W0                                          (whole-space power, direct field)
%            + omega/(8 pi^2) PV int int Im S[u_sc] k dk dphi   (regular part, k < max branch pt)
%            + sum_poles Wi.Wsrc                           (Sokhotski-Plemelj pole terms)
%   Beyond the last branch point Im S[u_sc] = 0 for real omega except at the poles.
[phi, wphi] = lamb.quad.azimuthNodes(opts.NPhi);
bp = iface.branchPoints(omega);
[k, wk] = lamb.quad.windowNodes([0, bp], opts.NGauss, lamb.quad.peaksFor(iface, omega, poles, leaky));
amp = lamb.incidentAmplitudes(src, iface, k, phi, omega);
x   = iface.solve(amp, k, omega);
S   = lamb.energy.sourceFunctional(src, iface, k, phi, omega, x);
Ws.W0      = W0;
Ws.regular = omega/(8*pi^2)*sum(wk.*k.*(imag(S)*wphi.'));
Ws.pole    = [Wi.Wsrc];
Ws.total   = Ws.W0 + Ws.regular + sum(Ws.pole);
end
