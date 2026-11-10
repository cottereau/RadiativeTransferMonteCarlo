function [Wb, W0, Wleak] = bodyWaves(iface, src, omega, poles, leaky, opts)
%BODYWAVES  Body-wave powers from the far-field vertical flux via Parseval (eq:Wbody3).
%   [Wb, W0, Wleak] = lamb.energy.bodyWaves(iface, src, omega, poles, leaky, opts)
%
%   W_m = (1/4pi^2) int_0^{2pi} dphi int_0^{k_m} E_m(k) |D_m(k,phi)|^2 k dk,
%   D_m = direct (away) + scattered, added coherently.
%   k: split at every branch point, graded around near-real leaky poles and toward branch
%   points that a real pole nearly touches (lamb.quad.peaksFor),
%   phi: exact uniform trapezoid (opts.NPhi).
%
%   Wb    : struct of powers [W], fields = names from lamb.energy.outgoingWaves
%   W0    : whole-space power of the direct field (up + down), for the balance check
%   Wleak : 1-by-numel(leaky); part of the body power within |k - a| < min(1e4*b, 1e-2*a)
%           of each near-real leaky pole (the Lorentzian peak, quasi-trapped energy;
%           the window misses 0.006 % of the peak and is approximate when b is not << a)
[phi, wphi] = lamb.quad.azimuthNodes(opts.NPhi);
bp    = iface.branchPoints(omega);
peaks = lamb.quad.peaksFor(iface, omega, poles, leaky);
waves = lamb.energy.outgoingWaves(iface, src.zs >= 0);
W0 = 0;  Wleak = zeros(1, numel(leaky));
for w = waves
    km = w.med.wavenumber(omega, w.type);
    [k, wk] = lamb.quad.windowNodes([0, bp(bp < km), km], opts.NGauss, peaks);
    amp = lamb.incidentAmplitudes(src, iface, k, phi, omega, true);
    x   = iface.solve(amp, k, omega);
    D   = x.(w.field);
    if ~isempty(w.away), D = D + amp.(w.away); end
    g   = wk.*k.*lamb.energy.fluxDensity(w.med, k, omega, w.type)/(4*pi^2);   % N-by-1
    dW  = g.*(abs(D).^2*wphi.');                                            % per node
    Wb.(w.name) = sum(dW);
    for j = 1:numel(leaky)
        Wleak(j) = Wleak(j) + sum(dW(abs(k - leaky(j).a) < min(1e4*leaky(j).b, 1e-2*leaky(j).a)));
    end
    if ~isempty(w.away)
        W0 = W0 + sum(g.*((abs(amp.(w.away)).^2 + abs(amp.(w.toward)).^2)*wphi.'));
    end
end
end
