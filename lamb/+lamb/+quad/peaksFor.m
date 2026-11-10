function peaks = peaksFor(iface, omega, poles, leaky)
%PEAKSFOR  Near-singular points of the real-axis integrands, as [a b] rows for windowNodes.
%   peaks = lamb.quad.peaksFor(iface, omega, poles, leaky)
%
%   - near-real leaky pole k_L = a + 1i*b inside a window: a Lorentzian of half-width b.
%   - real pole k_p within 1 % above the largest branch point k_n (weakly bound mode, e.g.
%     Stoneley with c -> min(beta)): the integrand near k_n varies on the scale
%     delta = k_p - k_n, so nodes are graded toward k_n with [a b] = [k_n, delta]
%
kn = max(iface.branchPoints(omega));
peaks = zeros(0, 2);
for j = 1:numel(leaky)
    peaks(end+1, :) = [leaky(j).a, leaky(j).b];                 %#ok<AGROW>
end
for j = 1:numel(poles)
    delta = poles(j).k - kn;
    if delta < 1e-2*kn
        peaks(end+1, :) = [kn, delta];                          %#ok<AGROW>
    end
end
end
