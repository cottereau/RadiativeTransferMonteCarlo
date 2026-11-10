function [r, lT, den] = poleNullVectors(iface, kp, omega)
%POLENULLVECTORS  Right/left null vectors of M(k_p) and the residue denominator (eq:res).
%   [r, lT, den] = lamb.poleNullVectors(iface, kp, omega)
%   Res_{k_p} M^{-1} = r*lT/den,  den = lT*M'(k_p)*r,  r and lT' unit-norm (svd).
%
%   M' is a centred difference with Richardson extrapolation (error O(h^4)). The step is well
%   below the distance delta from k_p to the nearest branch point, because nu ~ sqrt(k - k_n)
%   varies on that scale: a plain centred difference has error ~ (h/delta)^2/6, which is
%   1.7e-7 for a weakly bound Stoneley mode (delta = 1e-4 k_p) with h = 1e-3 delta.
M = iface.matrix(kp, omega);
[U, ~, V] = svd(M);
r  = V(:, end);
lT = U(:, end)';
h  = min(1e-4*kp, 1e-2*(kp - max(iface.branchPoints(omega))));
D  = @(h) (iface.matrix(kp + h, omega) - iface.matrix(kp - h, omega))/(2*h);
dM = (4*D(h/2) - D(h))/3;
den = lT*dM*r;
end
