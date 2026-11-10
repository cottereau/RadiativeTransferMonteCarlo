function M = stoneleyModel(vp1, vs1, rho1, vp2, vs2, rho2, freq, het1, het2)
%STONELEYMODEL  Analytical Stoneley interface wave model with scattering cross-sections.
%
%   M = stoneleyModel(vp1, vs1, rho1, vp2, vs2, rho2, freq, het1, het2)
%
%   Calculates Stoneley wave existence, phase velocity, displacement eigenfunctions,
%   energy partition, and complete scattering cross-sections between Solid 1 (z >= 0)
%   and Solid 2 (z <= 0) in welded contact.
%
%   Mathematical Foundation:
%     The Stoneley wave (Stoneley, 1924; Cagniard, 1962) is a true bound
%     interface wave between two elastic solids. It exists ONLY for a
%     restricted range of acoustic and density contrasts.
%
%     The boundary conditions at z = 0 are:
%       1. u_x^{(1)}(0) = u_x^{(2)}(0)       (continuity of tangential displacement)
%       2. u_z^{(1)}(0) = u_z^{(2)}(0)       (continuity of normal displacement)
%       3. tau_zz^{(1)}(0) = tau_zz^{(2)}(0) (continuity of normal traction)
%       4. tau_xz^{(1)}(0) = tau_xz^{(2)}(0) (continuity of shear traction)
%
%     This yields a 4x4 linear system M_mat(c) * [A1; A2; B1; B2] = 0.
%     Its determinant, normalized by the static (c/c_limit)^4 factor, is a
%     real function that is negative everywhere below the Stoneley root and
%     changes sign at most once in (0, min(vs1, vs2)) (Cagniard, 1962).  A
%     bound mode therefore exists if and only if that normalized determinant
%     is positive at the upper bound -- an exact sign test, with no residual
%     threshold and no lower cut-off.
%
%     Scattering Cross-Sections (Snieder 1986, Xu et al. 2021 extension):
%       Modal normalization:
%         Inorm = int_0^inf (|ux1|^2 + |uz1|^2) dz + (rho2/rho1) int_{-inf}^0 (|ux2|^2 + |uz2|^2) dz
%
%       In-plane 2D scattering (Stoneley -> Stoneley):
%         sig_StSt = w^3 / (8 * c_St^3 * Inorm^2) * [ Ir4_1 * e2_1 + (rho2/rho1)^2 * Ir4_2 * e2_2 ]
%
%       Bulk leakage radiation into Solid 1:
%         sig_StP1  = w*kP1^3 / (12*pi*c_St) * [Inorm_1 / Inorm] * e2_1
%         sig_StSv1 = w*kS1^3 / (24*pi*c_St) * [Ir14r2_1 / Inorm] * e2_1
%         sig_StSh1 = w*kS1^3 / ( 8*pi*c_St) * [Ir1_1 / Inorm] * e2_1
%
%       Bulk leakage radiation into Solid 2:
%         sig_StP2  = w*kP2^3 / (12*pi*c_St) * [(rho2/rho1)*Inorm_2 / Inorm] * e2_2
%         sig_StSv2 = w*kS2^3 / (24*pi*c_St) * [(rho2/rho1)*Ir14r2_2 / Inorm] * e2_2
%         sig_StSh2 = w*kS2^3 / ( 8*pi*c_St) * [(rho2/rho1)*Ir1_2 / Inorm] * e2_2
%
%       Total cross section:
%         Sigma_total = n1*(sig_StSt_1 + sig_StP1 + sig_StSv1 + sig_StSh1) + ...
%                       n2*(sig_StSt_2 + sig_StP2 + sig_StSv2 + sig_StSh2)
%         meanFreePath = 1 / Sigma_total
%
%   Inputs:
%     vp1, vs1, rho1 : Solid 1 parameters (z >= 0)
%     vp2, vs2, rho2 : Solid 2 parameters (z <= 0)
%     freq           : Frequency in Hz (default: 1000)
%     het1, het2     : Optional heterogeneity structs with fields:
%                      .V (volume, m^3), .drho (delta_rho/rho), .n (density, m^-3)
%
%   Outputs in struct M:
%     .exists        : Boolean flag (true if bound Stoneley mode exists)
%     .c_Stoneley    : Phase velocity (m/s), or NaN if non-existent
%     .searchLowerBound, .searchUpperBound : Velocity bracket used for the
%                         determinant sign test (m/s)
%     .boundaryResidual : Relative residual of the boundary matrix, at the
%                         root when it exists and at the upper bound otherwise
%     .k_Stoneley    : Wavenumber w/c (rad/m)
%     .lambda_St     : Wavelength (m)
%     .u_x1, .u_z1   : Displacement handles in solid 1 (z >= 0)
%     .u_x2, .u_z2   : Displacement handles in solid 2 (z <= 0)
%     .W1, .W2       : Energy density handles in solid 1 and solid 2
%     .E1, .E2       : Integrated energy in solid 1 and solid 2
%     .frac1         : Energy fraction in solid 1 (eta_1)
%     .frac2         : Energy fraction in solid 2 (eta_2)
%     .Inorm         : Modal energy normalization integral
%     .Inorm_1, .Inorm_2 : Per-medium parts of Inorm (medium 2 unweighted)
%     .U_Stoneley    : Group velocity (= phase velocity, non-dispersive here)
%     .c_P1, .c_S1, .c_P2, .c_S2 : bulk-to-Stoneley coefficients (m^4)
%     .sigBar_P1St, .sigBar_Sv1St, .sigBar_Sh1St,
%     .sigBar_P2St, .sigBar_Sv2St, .sigBar_Sh2St :
%                      angle-averaged, depth-integrated bulk-to-Stoneley
%                      cross sections (m^3); see the header note on units
%     .lame1, .lame2 : delta-lambda and delta-mu channels per medium, both
%                      directions, from interfaceLameChannels. Set
%                      het.dlam / het.dmu to activate; zero by default.
%     .sig_StSt      : In-plane Stoneley-to-Stoneley cross section (m^2)
%     .sig_StSt_1, .sig_StSt_2 : Per-medium parts of sig_StSt (m^2)
%     .sig_StP1      : Stoneley-to-P1 cross section (m^2)
%     .sig_StSv1     : Stoneley-to-Sv1 cross section (m^2)
%     .sig_StSh1     : Stoneley-to-Sh1 cross section (m^2)
%     .sig_StP2      : Stoneley-to-P2 cross section (m^2)
%     .sig_StSv2     : Stoneley-to-Sv2 cross section (m^2)
%     .sig_StSh2     : Stoneley-to-Sh2 cross section (m^2)
%     .Sigma_1, .Sigma_2 : Per-medium scattering rates (1/m)
%     .Sigma_total   : Total scattering rate (1/m)
%     .meanFreePath  : Mean free path ell = 1/Sigma_total (m)
%     .meanFreeTime  : Mean free time tau = ell / c_Stoneley (s)

if nargin < 7 || isempty(freq), freq = 1000; end
validateStoneleyInputs(vp1, vs1, rho1, vp2, vs2, rho2, freq);
w = 2 * pi * freq;

% Parse heterogeneity parameters (default: 5% fluctuations, unit scatterer)
if nargin < 8 || isempty(het1)
    het1 = struct('V', 1.0, 'drho', 0.05, 'n', 1.0);
elseif isnumeric(het1)
    het1 = struct('V', sqrt(het1), 'drho', 1.0, 'n', 1.0);
end
if nargin < 9 || isempty(het2)
    het2 = struct('V', 1.0, 'drho', 0.05, 'n', 1.0);
elseif isnumeric(het2)
    het2 = struct('V', sqrt(het2), 'drho', 1.0, 'n', 1.0);
end
het1 = withLameDefaults(het1);
het2 = withLameDefaults(het2);
validateHeterogeneity(het1, 'het1');
validateHeterogeneity(het2, 'het2');

e2_1 = het1.V^2 * het1.drho^2;
e2_2 = het2.V^2 * het2.drho^2;
n1   = het1.n;
n2   = het2.n;

c_limit = min(vs1, vs2);
% Keep clear of both endpoints: c = c_limit makes the shear decay rate vanish,
% and below ~1e-4*c_limit the 4x4 determinant loses the sign of its static
% branch to round-off.  Neither is a physical cut-off.
c_min = 1e-4 * c_limit;
c_max = c_limit - 8 * eps(c_limit);

residualFunc = @(c) computeStoneleyResidual(c, w, vp1, vs1, rho1, vp2, vs2, rho2);
detFunc = @(c) computeStoneleyDet(c, w, vp1, vs1, rho1, vp2, vs2, rho2, c_limit);

% Existence is the sign of the normalized determinant at the upper bound.
% The previous smallest-singular-value test could not do this: sigma_min
% decays like sqrt(c_limit-c) at the endpoint whether or not a mode exists,
% so a fixed tolerance was a coincidence -- two identical solids landed at
% 1.4e-8 against a 1e-8 threshold.
det_upper = detFunc(c_max);
if ~isfinite(det_upper) || det_upper <= 0
    M = emptyStoneleyResult();
    M.message = sprintf(['No bound Stoneley wave exists for this pair: the ', ...
        'normalized boundary determinant is %+.3e at %.6g m/s and does not ', ...
        'change sign above %.6g m/s.'], det_upper, c_max, c_min);
    M.boundaryResidual = residualFunc(c_max);
    M.searchLowerBound = c_min;
    M.searchUpperBound = c_max;
    return;
end

det_lower = detFunc(c_min);
if ~isfinite(det_lower) || det_lower >= 0
    error('stoneleyModel:UnbracketedRoot', ...
        ['The Stoneley determinant is %+.3e at %.6g m/s, so the root is not ', ...
         'bracketed by [%.6g, %.6g] m/s.'], det_lower, c_min, c_min, c_max);
end
c_Stoneley = fzero(detFunc, [c_min, c_max]);

% Compute the null vector from the scaled matrix.  Row scaling does not alter
% the null space and avoids traction rows dominating displacement rows.
[best_residual, M_mat, M_scaled] = computeStoneleyResidual(c_Stoneley, w, vp1, vs1, rho1, vp2, vs2, rho2);
[~, ~, V_svd] = svd(M_scaled);
coeffs = V_svd(:, 4); % Null vector
A1 = coeffs(1);
A2 = coeffs(2);
B1 = coeffs(3);
B2 = coeffs(4);

k = w / c_Stoneley;
bp1 = sqrt(k^2 - (w/vp1)^2);
bs1 = sqrt(k^2 - (w/vs1)^2);
bp2 = sqrt(k^2 - (w/vp2)^2);
bs2 = sqrt(k^2 - (w/vs2)^2);

% Displacements:
% Solid 1 (z >= 0):
ux1 = @(z) (1i*k*A1*exp(-bp1*z) + bs1*A2*exp(-bs1*z));
uz1 = @(z) (-bp1*A1*exp(-bp1*z) + 1i*k*A2*exp(-bs1*z));

% Solid 2 (z <= 0):
ux2 = @(z) (1i*k*B1*exp(bp2*z) - bs2*B2*exp(bs2*z));
uz2 = @(z) (bp2*B1*exp(bp2*z) + 1i*k*B2*exp(bs2*z));

% Time-averaged elastic energy densities: kinetic plus strain energy.
dux1_dz = @(z) (-1i*k*bp1*A1*exp(-bp1*z) - bs1^2*A2*exp(-bs1*z));
duz1_dz = @(z) ( bp1^2*A1*exp(-bp1*z) - 1i*k*bs1*A2*exp(-bs1*z));
dux2_dz = @(z) ( 1i*k*bp2*B1*exp( bp2*z) - bs2^2*B2*exp( bs2*z));
duz2_dz = @(z) ( bp2^2*B1*exp( bp2*z) + 1i*k*bs2*B2*exp( bs2*z));
W1 = elasticEnergyDensity(ux1, uz1, dux1_dz, duz1_dz, k, rho1, vp1, vs1, w);
W2 = elasticEnergyDensity(ux2, uz2, dux2_dz, duz2_dz, k, rho2, vp2, vs2, w);

% Every eigenfunction is a sum of decaying exponentials, so the semi-infinite
% depth integrals converge and integral() transforms them itself.  A fixed
% 30/min(b) window is both more code and less accurate: near the existence
% threshold it spans ~1e8 m and the adaptive quadrature then drifts 2.6e-9
% from the exact closed form, against 2e-14 here.
E1 = integral(W1, 0, Inf, 'RelTol', 1e-10, 'AbsTol', 0);
E2 = integral(W2, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);

E_tot = E1 + E2;
frac1 = E1 / E_tot;
frac2 = E2 / E_tot;

% ---- Exact Modal Depth Integrals for Scattering Cross-Sections ----
r1_1 = @(z) abs(ux1(z));
r2_1 = @(z) abs(uz1(z));
r1_2 = @(z) abs(ux2(z));
r2_2 = @(z) abs(uz2(z));

Inorm_1 = integral(@(z) r1_1(z).^2 + r2_1(z).^2, 0, Inf, 'RelTol', 1e-10, 'AbsTol', 0);
Inorm_2 = integral(@(z) r1_2(z).^2 + r2_2(z).^2, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);
Inorm   = Inorm_1 + (rho2 / rho1) * Inorm_2;

% Quartic integrals
Ir4_1 = integral(@(z) r1_1(z).^4 + 2*r2_1(z).^4, 0, Inf, 'RelTol', 1e-10, 'AbsTol', 0);
Ir4_2 = integral(@(z) r1_2(z).^4 + 2*r2_2(z).^4, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);

% Quadratic components for bulk radiation
Ir1_1    = integral(@(z) r1_1(z).^2, 0, Inf, 'RelTol', 1e-10, 'AbsTol', 0);
Ir2_1    = integral(@(z) r2_1(z).^2, 0, Inf, 'RelTol', 1e-10, 'AbsTol', 0);
Ir14r2_1 = Ir1_1 + 4 * Ir2_1;

Ir1_2    = integral(@(z) r1_2(z).^2, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);
Ir2_2    = integral(@(z) r2_2(z).^2, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);
Ir14r2_2 = Ir1_2 + 4 * Ir2_2;

kP1 = w / vp1;  kS1 = w / vs1;
kP2 = w / vp2;  kS2 = w / vs2;

% 1. In-plane Stoneley -> Stoneley
sig_StSt_1 = (w^3 / (8 * c_Stoneley^3 * Inorm^2)) * Ir4_1 * e2_1;
sig_StSt_2 = (w^3 / (8 * c_Stoneley^3 * Inorm^2)) * ((rho2/rho1)^2 * Ir4_2) * e2_2;
sig_StSt   = sig_StSt_1 + sig_StSt_2;

% 2. Bulk radiation into Solid 1 (P, Sv, Sh)
% The P kernel is (r1^2 + r2^2) (Xu et al., 2021, B.6), so its depth integral
% is the per-medium part of Inorm.  A single halfspace makes that ratio 1,
% which is why B.17 shows none; here Inorm_1/Inorm and (rho2/rho1)*Inorm_2/Inorm
% share one normalization and must sum to 1.
sig_StP1  = (w * kP1^3 / (12 * pi * c_Stoneley)) * (Inorm_1 / Inorm) * e2_1;
sig_StSv1 = (w * kS1^3 / (24 * pi * c_Stoneley)) * (Ir14r2_1 / Inorm) * e2_1;
sig_StSh1 = (w * kS1^3 / ( 8 * pi * c_Stoneley)) * (Ir1_1 / Inorm) * e2_1;

% 3. Bulk radiation into Solid 2 (P, Sv, Sh)
sig_StP2  = (w * kP2^3 / (12 * pi * c_Stoneley)) * ((rho2/rho1) * Inorm_2 / Inorm) * e2_2;
sig_StSv2 = (w * kS2^3 / (24 * pi * c_Stoneley)) * ((rho2/rho1) * Ir14r2_2 / Inorm) * e2_2;
sig_StSh2 = (w * kS2^3 / ( 8 * pi * c_Stoneley)) * ((rho2/rho1) * Ir1_2 / Inorm) * e2_2;

% ---- Inverse channels: bulk -> Stoneley ------------------------------
% Reverse of the six bulk channels above, i.e. the analogue of
% RayleighModel's cP/cS and Xu et al. B.3-B.5. For a bulk wave of mode b in
% medium j arriving at incidence psi from the interface normal, scattering
% off a density heterogeneity at depth z:
%
%   sigma_{P ->St}(z,psi) = c_Pj*[ r1_j(z)^2 sin^2(psi) + 2 r2_j(z)^2 cos^2(psi) ]
%   sigma_{Sv->St}(z,psi) = c_Sj*[ r1_j(z)^2 cos^2(psi) + 2 r2_j(z)^2 sin^2(psi) ]
%   sigma_{Sh->St}(z)     = c_Sj*  r1_j(z)^2
%
% Units: c_* are m^4, kernels m^-2, so each sigma is m^2 AT A GIVEN DEPTH --
% unlike the guided->bulk channels above, which are already integrated over
% the scatterer depth. The two directions are neither equal nor
% interchangeable; what relates them is the density of states, asserted as a
% detailed-balance check in testInterfaceWavesAndWedge.m.
%
% The (rho_2/rho_1) weight is the same one the forward channels carry. It has
% to be: it cancels in detailed balance together with Inorm_j, leaving the
% medium-independent ratio pi*v^3/(w*c*U).
U_St = c_Stoneley;                      % no dispersion between two half-spaces
c_P1 = (w^3 / (8 * vp1 * c_Stoneley * U_St * Inorm)) * e2_1;
c_S1 = (w^3 / (8 * vs1 * c_Stoneley * U_St * Inorm)) * e2_1;
c_P2 = (w^3 / (8 * vp2 * c_Stoneley * U_St * Inorm)) * (rho2/rho1) * e2_2;
c_S2 = (w^3 / (8 * vs2 * c_Stoneley * U_St * Inorm)) * (rho2/rho1) * e2_2;

% Angle-averaged and depth-integrated (the B.6-B.8 analogue).
sigBar_P1St  = (2/3) * c_P1 * Inorm_1;
sigBar_Sv1St = (1/3) * c_S1 * (Ir1_1 + 4*Ir2_1);
sigBar_Sh1St =         c_S1 * Ir1_1;
sigBar_P2St  = (2/3) * c_P2 * Inorm_2;
sigBar_Sv2St = (1/3) * c_S2 * (Ir1_2 + 4*Ir2_2);
sigBar_Sh2St =         c_S2 * Ir1_2;

% ---- delta-lambda and delta-mu channels -------------------------------
% Both directions again, for heterogeneity in the Lame parameters of each
% solid. Zero unless het.dlam / het.dmu were set.
K1_lame = interfaceLameKernels(ux1, uz1, dux1_dz, duz1_dz, k);
K2_lame = interfaceLameKernels(ux2, uz2, dux2_dz, duz2_dz, k);
C1_lame = interfaceLameChannels(K1_lame, [0 Inf], w, vp1, vs1, rho1, rho1, ...
    c_Stoneley, U_St, Inorm, het1.V, het1.dlam, het1.dmu);
C2_lame = interfaceLameChannels(K2_lame, [-Inf 0], w, vp2, vs2, rho2, rho1, ...
    c_Stoneley, U_St, Inorm, het2.V, het2.dlam, het2.dmu);

% Total cross section Sigma (1/m)
Sigma_1 = n1 * (sig_StSt_1 + sig_StP1 + sig_StSv1 + sig_StSh1 ...
    + C1_lame.sig_lam_P + C1_lame.sig_mu_P + C1_lame.sig_mu_Sv + C1_lame.sig_mu_Sh);
Sigma_2 = n2 * (sig_StSt_2 + sig_StP2 + sig_StSv2 + sig_StSh2 ...
    + C2_lame.sig_lam_P + C2_lame.sig_mu_P + C2_lame.sig_mu_Sv + C2_lame.sig_mu_Sh);
Sigma_total = Sigma_1 + Sigma_2;

meanFreePath = 1.0 / Sigma_total;
meanFreeTime = meanFreePath / c_Stoneley;

% Effective scattering proxy:
eff_sigma = @(s1, s2) frac1 * s1 + frac2 * s2;

% Populate struct
M.exists        = true;
M.c_Stoneley    = c_Stoneley;
M.boundaryResidual = best_residual;
M.searchLowerBound = c_min;
M.searchUpperBound = c_max;
M.k_Stoneley    = k;
M.lambda_St     = 2*pi / k;
M.bp1           = bp1;
M.bs1           = bs1;
M.bp2           = bp2;
M.bs2           = bs2;
M.u_x1          = ux1;
M.u_z1          = uz1;
M.u_x2          = ux2;
M.u_z2          = uz2;
M.W1            = W1;
M.W2            = W2;
M.E1            = E1;
M.E2            = E2;
M.frac1         = frac1;
M.frac2         = frac2;
M.Inorm         = Inorm;
M.Inorm_1       = Inorm_1;
M.Inorm_2       = Inorm_2;
M.U_Stoneley    = U_St;
M.c_P1          = c_P1;   M.c_S1 = c_S1;
M.c_P2          = c_P2;   M.c_S2 = c_S2;
M.sigBar_P1St   = sigBar_P1St;
M.sigBar_Sv1St  = sigBar_Sv1St;
M.sigBar_Sh1St  = sigBar_Sh1St;
M.sigBar_P2St   = sigBar_P2St;
M.sigBar_Sv2St  = sigBar_Sv2St;
M.sigBar_Sh2St  = sigBar_Sh2St;
M.lame1         = C1_lame;
M.lame2         = C2_lame;
M.sig_StSt      = sig_StSt;
M.sig_StSt_1    = sig_StSt_1;
M.sig_StSt_2    = sig_StSt_2;
M.sig_StP1      = sig_StP1;
M.sig_StSv1     = sig_StSv1;
M.sig_StSh1     = sig_StSh1;
M.sig_StP2      = sig_StP2;
M.sig_StSv2     = sig_StSv2;
M.sig_StSh2     = sig_StSh2;
M.Sigma_1       = Sigma_1;
M.Sigma_2       = Sigma_2;
M.Sigma_total   = Sigma_total;
M.meanFreePath  = meanFreePath;
M.meanFreeTime  = meanFreeTime;
M.eff_sigma     = eff_sigma;

end

function [det_val, M] = computeStoneleyDet(c, w, vp1, vs1, rho1, vp2, vs2, rho2, c_limit)
    k = w / c;
    bp1 = sqrt(k^2 - (w/vp1)^2);
    bs1 = sqrt(k^2 - (w/vs1)^2);
    bp2 = sqrt(k^2 - (w/vp2)^2);
    bs2 = sqrt(k^2 - (w/vs2)^2);
    
    mu1 = rho1 * vs1^2;
    mu2 = rho2 * vs2^2;
    ks1_sq = (w/vs1)^2;
    ks2_sq = (w/vs2)^2;
    
    M = zeros(4, 4);
    
    % Row 1: u_x continuity (ux1 - ux2 = 0)
    M(1, 1) = 1i * k;
    M(1, 2) = bs1;
    M(1, 3) = -1i * k;
    M(1, 4) = bs2;
    
    % Row 2: u_z continuity (uz1 - uz2 = 0)
    M(2, 1) = -bp1;
    M(2, 2) = 1i * k;
    M(2, 3) = -bp2;
    M(2, 4) = -1i * k;
    
    % Row 3: tau_zz continuity (tau_zz1 - tau_zz2 = 0)
    M(3, 1) = mu1 * (2*k^2 - ks1_sq);
    M(3, 2) = -2i * mu1 * k * bs1;
    M(3, 3) = -mu2 * (2*k^2 - ks2_sq);
    M(3, 4) = -2i * mu2 * k * bs2;
    
    % Row 4: tau_xz continuity (tau_xz1 - tau_xz2 = 0)
    M(4, 1) = -2i * mu1 * k * bp1;
    M(4, 2) = -mu1 * (2*k^2 - ks1_sq);
    M(4, 3) = -2i * mu2 * k * bp2;
    M(4, 4) = mu2 * (2*k^2 - ks2_sq);
    
    % Normalized determinant.  det(M) is exactly real for this matrix.  It
    % vanishes as (c/c_limit)^4 at c -> 0; dividing that static factor out
    % leaves an O(1) function whose only zero is the Stoneley root, the same
    % device computeSecular uses for the trivial c = 0 root of Rayleigh.
    scale = (mu1 + mu2)^2 * (k^6);
    if nargin >= 9
        scale = scale * (c / c_limit)^4;
    end
    det_val = real(det(M)) / scale;
end

function [residual, M, M_scaled] = computeStoneleyResidual(c, w, vp1, vs1, rho1, vp2, vs2, rho2)
    [~, M] = computeStoneleyDet(c, w, vp1, vs1, rho1, vp2, vs2, rho2, 1);
    k = w / c;
    mu1 = rho1 * vs1^2;
    mu2 = rho2 * vs2^2;

    % Rows 1--2 express displacement and rows 3--4 traction.  Scaling each
    % family puts them on comparable dimensionless magnitudes before SVD.
    row_scale = [k; k; (mu1 + mu2) * k^2; (mu1 + mu2) * k^2];
    M_scaled = diag(1 ./ row_scale) * M;
    singular_values = svd(M_scaled);
    residual = singular_values(end) / singular_values(1);
end

function M = emptyStoneleyResult()
    % Keep the public struct schema stable when the bound mode is absent.
    M = struct( ...
        'exists', false, ...
        'c_Stoneley', NaN, ...
        'searchLowerBound', NaN, ...
        'searchUpperBound', NaN, ...
        'k_Stoneley', NaN, ...
        'lambda_St', NaN, ...
        'bp1', NaN, ...
        'bs1', NaN, ...
        'bp2', NaN, ...
        'bs2', NaN, ...
        'u_x1', [], ...
        'u_z1', [], ...
        'u_x2', [], ...
        'u_z2', [], ...
        'W1', [], ...
        'W2', [], ...
        'E1', NaN, ...
        'E2', NaN, ...
        'frac1', NaN, ...
        'frac2', NaN, ...
        'Inorm', NaN, ...
        'Inorm_1', NaN, ...
        'Inorm_2', NaN, ...
        'U_Stoneley', NaN, ...
        'c_P1', NaN, 'c_S1', NaN, 'c_P2', NaN, 'c_S2', NaN, ...
        'sigBar_P1St', NaN, 'sigBar_Sv1St', NaN, 'sigBar_Sh1St', NaN, ...
        'sigBar_P2St', NaN, 'sigBar_Sv2St', NaN, 'sigBar_Sh2St', NaN, ...
        'lame1', [], 'lame2', [], ...
        'sig_StSt', NaN, ...
        'sig_StSt_1', NaN, ...
        'sig_StSt_2', NaN, ...
        'sig_StP1', NaN, ...
        'sig_StSv1', NaN, ...
        'sig_StSh1', NaN, ...
        'sig_StP2', NaN, ...
        'sig_StSv2', NaN, ...
        'sig_StSh2', NaN, ...
        'Sigma_1', NaN, ...
        'Sigma_2', NaN, ...
        'Sigma_total', NaN, ...
        'meanFreePath', NaN, ...
        'meanFreeTime', NaN, ...
        'eff_sigma', [], ...
        'boundaryResidual', NaN, ...
        'message', '');
end

function W = elasticEnergyDensity(ux, uz, dux_dz, duz_dz, k, rho, vp, vs, w)
    % Time-averaged energy density for a 2-D plane-strain elastic phasor.
    mu = rho * vs^2;
    lambda = rho * vp^2 - 2 * mu;
    exx = @(z) 1i * k * ux(z);
    ezz = @(z) duz_dz(z);
    exz = @(z) 0.5 * (1i * k * uz(z) + dux_dz(z));
    K = @(z) 0.25 * rho * w^2 * (abs(ux(z)).^2 + abs(uz(z)).^2);
    U = @(z) 0.25 * lambda * abs(exx(z) + ezz(z)).^2 + ...
             0.5 * mu * (abs(exx(z)).^2 + abs(ezz(z)).^2 + 2 * abs(exz(z)).^2);
    W = @(z) K(z) + U(z);
end

function validateStoneleyInputs(vp1, vs1, rho1, vp2, vs2, rho2, freq)
    validatePositiveScalar(vp1, 'vp1');
    validatePositiveScalar(vs1, 'vs1');
    validatePositiveScalar(rho1, 'rho1');
    validatePositiveScalar(vp2, 'vp2');
    validatePositiveScalar(vs2, 'vs2');
    validatePositiveScalar(rho2, 'rho2');
    validatePositiveScalar(freq, 'freq');
    if vp1 <= vs1 || vp2 <= vs2
        error('stoneleyModel:InvalidInput', ...
            'Each isotropic elastic solid must satisfy vp > vs.');
    end
end

function het = withLameDefaults(het)
    % Optional, default zero: every existing caller keeps its cross sections.
    if ~isfield(het, 'dlam') || isempty(het.dlam), het.dlam = 0; end
    if ~isfield(het, 'dmu')  || isempty(het.dmu),  het.dmu  = 0; end
end

function validateHeterogeneity(het, name)
    if ~isstruct(het) || ~isscalar(het) || ~all(isfield(het, {'V', 'drho', 'n'}))
        error('stoneleyModel:InvalidHeterogeneity', ...
            '%s must be a scalar struct with fields V, drho, and n.', name);
    end
    validateNonnegativeScalar(het.V, [name '.V']);
    validateFiniteRealScalar(het.drho, [name '.drho']);
    validateNonnegativeScalar(het.n, [name '.n']);
    validateFiniteRealScalar(het.dlam, [name '.dlam']);
    validateFiniteRealScalar(het.dmu,  [name '.dmu']);
end

function validatePositiveScalar(value, name)
    validateFiniteRealScalar(value, name);
    if value <= 0
        error('stoneleyModel:InvalidInput', '%s must be a positive finite real scalar.', name);
    end
end

function validateNonnegativeScalar(value, name)
    validateFiniteRealScalar(value, name);
    if value < 0
        error('stoneleyModel:InvalidInput', '%s must be a non-negative finite real scalar.', name);
    end
end

function validateFiniteRealScalar(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('stoneleyModel:InvalidInput', '%s must be a finite real scalar.', name);
    end
end
