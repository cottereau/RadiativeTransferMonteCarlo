function [sigma, vEff, alphaEff] = calcSigma(mat)
%% calcSigma - DSCS from EXACT scattering by a single sphere (partial-wave T-matrix)
%
%   [sigma, vEff, alphaEff] = TMatrix.calcSigma(mat)
%
% Used by MaterialClass.CalcSigma when mat.SpectralLaw = 'tmatrix'. Valid for
% spherical inclusions only. Replaces the Born approximation of
% calcSigmaIsotropicCorrelation in the strong-contrast / resonant regime.
%
% Motivation: Born (Ryzhik et al. 1996, doi:10.1016/S0165-2125(96)00021-2)
% assumes weak fluctuations. In the
% two-phase case of the Manda data set the coefficients of variation of
% lambda and mu are 42% and 50%, and k*Lc ranges from 1.3 to 3.5 (resonant
% region). Measured: Born overestimates the frequency dependence of Sigma*
% by a factor ~4 (+69% from 300 to 800 kHz, against +17% for SEM).
%
% For each pair of modes,
%
%   sigma_ij(theta) = v_i * n0 * |f_ij(theta)|^2 * S(q_ij(theta))
%
% The factor v_i (velocity of the incident mode) converts from 1/m to 1/s,
% the convention used for sigma in MaterialClass (prepareSigma sets
% meanFreeTime = 1/Sigma).
%
% q depends on the mode pair:
%   q_PP = 2 kP sin(th/2)   q_SS = 2 kS sin(th/2)
%   q_PS = sqrt(kP^2 + kS^2 - 2 kP kS cos th)
%
% Required fields of mat.SpectralParam:
%   SphereRadius          sphere radius a [m]
%   SphereVolumeFraction  volume fraction of the spheres
%   HostPhase             struct (lam, mu, rho) of the host matrix
%   InclusionPhase        struct (lam, mu, rho) of the sphere
% Optional:
%   StructureFactor       handle S(q); empty/absent = 1 (independent
%                         scatterers), which underestimates the suppression
%                         of forward scattering at high volume fraction.
%
% Outputs:
%   sigma     2x2 cell of handles sigma{i,j}(theta) [1/s]
%   vEff      [vp vs] effective phase velocities (optical theorem), INFORMATIVE
%   alphaEff  [alpha_p alpha_s] scattering attenuation [1/m], INFORMATIVE

if mat.acoustics
    error('TMatrix.calcSigma is implemented only for the elastic case');
end
prm = mat.SpectralParam;
required = {'SphereRadius','SphereVolumeFraction','HostPhase','InclusionPhase'};
missing = required(~isfield(prm, required));
if isempty(missing)
    missing = required(cellfun(@(c) isempty(prm.(c)), required));
end
if ~isempty(missing)
    error('TMatrix.calcSigma needs SpectralParam fields: %s', strjoin(missing,', '));
end
if isempty(mat.Frequency)
    error('Set mat.Frequency before computing the scattering cross-section');
end

a   = prm.SphereRadius;
n0  = prm.SphereVolumeFraction/((4/3)*pi*a^3);   % number density [1/m^3]
om  = 2*pi*mat.Frequency(1);

% velocities of the HOST matrix (not the effective ones of the object):
% the scattering problem is a sphere embedded in the matrix
hp  = prm.HostPhase;
vpH = sqrt((hp.lam+2*hp.mu)/hp.rho);
vsH = sqrt(hp.mu/hp.rho);
kp  = om/vpH;  ks = om/vsH;

sc = TMatrix.sphereScattering(a, prm.HostPhase, prm.InclusionPhase, om);
th = sc.th;

if ~isfield(prm,'StructureFactor') || isempty(prm.StructureFactor)
    Sf = @(q) ones(size(q));
else
    Sf = prm.StructureFactor;
end

qPP = 2*kp*sin(th/2);
qSS = 2*ks*sin(th/2);
qPS = sqrt(kp^2 + ks^2 - 2*kp*ks*cos(th));

g = {vpH*n0*sc.dsig_PP.*Sf(qPP), vpH*n0*sc.dsig_PS.*Sf(qPS); ...
     vsH*n0*sc.dsig_SP.*Sf(qPS), vsH*n0*sc.dsig_SS.*Sf(qSS)};

% --- effective velocity from the optical theorem (INFORMATIVE ONLY) ---
% Foldy's formula for independent scatterers (first order in n0):
%
%   K^2 = k^2 + 4*pi*n0*f(0)
%
% with f the far-field amplitude normalised as u_sc ~ f(theta) exp(ikr)/r
% for a unit-amplitude incident plane wave (the normalisation returned by
% TMatrix.sphereScattering, see farFat there).
%
% Reference: L. L. Foldy, "The multiple scattering of waves. I. General
% theory of isotropic scattering by randomly distributed scatterers",
% Phys. Rev. 67 (1945) 107-119, doi:10.1103/PhysRev.67.107. In the notation of C. M. Linton and
% P. A. Martin, "Multiple scattering by multiple spheres: a new proof of
% the Lloyd-Berry formula for the effective wavenumber", SIAM J. Appl.
% Math. 66(5) (2006) 1649-1668, doi:10.1137/050636401, Foldy's formula is their Eq. (1.1); the
% same first-order term written with the far-field pattern is
% delta_1 = -(4*pi*i/k) f_LM(0) in their Eqs. (1.2), (1.4)-(1.5), with
% f_LM defined by their Eq. (1.3), u_sc ~ (ikr)^-1 exp(ikr) f_LM(Theta).
% Since f_LM = i*k*f, delta_1 = 4*pi*f(0), the expression above.
% The scalar formula is applied separately to the P and S coherent waves
% (forward amplitudes f_PP(0) and f_SS(0)); for a sphere the converted
% amplitudes f_PS and f_SP vanish at theta = 0.
%
% Re(K) gives the effective phase velocity, Im(K) the scattering
% attenuation. NOT used anywhere: transport keeps mat.vp and mat.vs as set
% by the user. It shows how far the effective velocities are from those of
% the matrix, the known inconsistency of this path (cross-section in the
% host matrix, transport at the effective velocity).
kEffP = sqrt(kp^2 + 4*pi*n0*sc.f_PP_fwd);
kEffS = sqrt(ks^2 + 4*pi*n0*sc.f_SS_fwd);
vEff     = [om/real(kEffP), om/real(kEffS)];
alphaEff = [imag(kEffP), imag(kEffS)];
fprintf(['TMatrix.calcSigma @ %.0f kHz | v_eff (optical theorem) = ' ...
         '%.0f / %.0f m/s  (matrix: %.0f / %.0f | object: %.0f / %.0f)\n' ...
         '                          scattering attenuation = ' ...
         '%.2f / %.2f 1/m\n'], ...
        mat.Frequency(1)/1e3, vEff, vpH, vsH, mat.vp, mat.vs, alphaEff);

% prepareSigmaOne integrates with adaptive quadrature, so sigma must be a
% handle that accepts vector theta
sigma = cell(2);
for i = 1:2
    for j = 1:2
        F = griddedInterpolant(th, g{i,j}, 'pchip', 'nearest');
        sigma{i,j} = @(x) F(x);
    end
end
end
