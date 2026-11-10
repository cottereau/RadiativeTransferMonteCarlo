function F = interfaceFormFactors(kI, kb, Phi, Lc, Ir1, Ir2, nDir)
%INTERFACEFORMFACTORS  Finite-size correction to the interface-mode channels.
%
%   F = interfaceFormFactors(kI, kb, Phi, Lc, Ir1, Ir2) returns the factor by
%   which a finite-size random heterogeneity changes each guided<->bulk
%   channel relative to a point scatterer.
%
%   The point-scatterer cross sections in scholteModel and stoneleyModel
%   carry V^2 (dm/m)^2, a constant. For a continuous random medium that
%   constant becomes the power spectral density of the fluctuation evaluated
%   at the momentum transfer,
%
%       q = k_b n_hat - k_I x_hat        (guided mode taken along +x)
%       |q|^2 = k_b^2 + k_I^2 - 2 k_b k_I sin(theta) cos(phi)
%
%   which no longer factors out of the angular integral: it turns an
%   isotropic point-scatterer pattern into an anisotropic phase function.
%   So each channel picks up
%
%       F = < kernel(theta,phi) * Phi(|q| Lc) > / < kernel(theta,phi) >
%
%   with the same angular kernels the cross sections were integrated from
%   (Xu et al. B.3-B.5 / B.9-B.11):
%
%       P  : r1^2 sin^2(theta) cos^2(phi) + r2^2 cos^2(theta)
%       Sv : r1^2 cos^2(theta) cos^2(phi) + r2^2 sin^2(theta)
%       Sh : r1^2 sin^2(phi)
%
%   and the depth integrals Ir1 = INT r1^2 dz, Ir2 = INT r2^2 dz supplying
%   the two weights, since the depth and angle parts separate.
%
%   Phi takes the DIMENSIONLESS argument q*Lc, matching MaterialClass's
%   convention (see the sigmaPP etc. handles there). Passing Phi = [] or a
%   constant returns F = 1 exactly, recovering the point-scatterer model.
%
%   WHAT THIS DOES NOT DO. It corrects the ANGULAR SHAPE only. Going fully
%   continuous also replaces the amplitude n*V^2*(dm/m)^2 by the PSD's own
%   normalisation; that substitution is a separate modelling choice and is
%   not made here. The evanescent vertical wavenumber of the guided mode is
%   likewise dropped from q -- valid while the decay rate is small against
%   |q|, which is the usual local-plane-wave assumption.
%
%   See also SCHOLTEMODEL, STONELEYMODEL, MATERIALCLASS/PREPARESIGMA.
%
%   References:
%     - Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Monte Carlo simulations of coupled
%       body- and Rayleigh-wave multiple scattering in elastic media, Geophys. J. Int. 228,
%       1213-1236. doi:10.1093/gji/ggab370

if nargin < 7 || isempty(nDir), nDir = 200000; end
if isempty(Phi) || isempty(Lc) || Lc == 0
    F = struct('P',1,'Sv',1,'Sh',1);
    return
end

% quasi-uniform directions on the sphere (Fibonacci)
ii = (0:nDir-1).' + 0.5;
ct = 1 - 2*ii/nDir;                       % cos(theta)
st = sqrt(max(0,1-ct.^2));
ga = pi*(1+sqrt(5));
ph = mod(ga*ii, 2*pi);

% momentum transfer between the guided mode (along +x) and the bulk wave
q = sqrt(max(0, kb^2 + kI^2 - 2*kb*kI*st.*cos(ph)));
w = Phi(q*Lc);

kerP  = Ir1*st.^2.*cos(ph).^2 + Ir2*ct.^2;
kerSv = Ir1*ct.^2.*cos(ph).^2 + Ir2*st.^2;
kerSh = Ir1*sin(ph).^2;

F.P  = mean(kerP .*w)/mean(kerP);
F.Sv = mean(kerSv.*w)/mean(kerSv);
F.Sh = mean(kerSh.*w)/mean(kerSh);
end
