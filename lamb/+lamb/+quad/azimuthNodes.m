function [phi, w] = azimuthNodes(N)
%AZIMUTHNODES  Uniform azimuths phi_j = 2*pi*j/N with weights 2*pi/N (row vectors).
%   The trapezoidal rule is exact for trigonometric polynomials of degree < N. Energy
%   integrands have degree <= 4 (moment tensor), so N >= 5 is exact.
phi = 2*pi*(0:N-1)/N;
w   = (2*pi/N)*ones(1, N);
end
