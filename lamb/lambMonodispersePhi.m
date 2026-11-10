function [Phi, Lc] = lambMonodispersePhi(eta, D)
%LAMBMONODISPERSEPHI  Normalised 3D Percus-Yevick PSDF of randomly placed hard spheres.
%   [Phi, Lc] = lambMonodispersePhi(eta, D)
%
%   eta  volume fraction,  D  sphere diameter (any length unit; Lc comes out in the same unit).
%   Thin wrapper: the spectral-route computation lives in MaterialClass.MonoDisperseSphere (3D).
m = MaterialClass;  m.d = 3;
Phi = m.MonoDisperseSphere(eta, D);
Lc = m.CorrelationLength;
end
