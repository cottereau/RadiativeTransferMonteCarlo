function M = scholteModel(vp_s, vs_s, rho_s, vf, rho_f, freq, het_s, het_f)
%SCHOLTEMODEL  Analytical Scholte interface wave model with scattering cross-sections.
%
%   M = scholteModel(vp_s, vs_s, rho_s, vf, rho_f, freq, het_s, het_f)
%
%   Calculates the exact Scholte phase velocity, displacement eigenfunctions,
%   energy density profile, solid/fluid energy partition, and the complete
%   modal scattering cross-sections (in-plane 2D scattering + 3D bulk leakage
%   into fluid and solid) for an elastic solid (z <= 0) and acoustic fluid (z >= 0).
%
%   Mathematical Foundation:
%     The Scholte wave (Scholte, 1947) is an interface acoustic-elastic wave
%     confined to the solid-fluid boundary. Its phase velocity c_Scholte is
%     the real root of the Scholte secular equation:
%
%       F(c) = 4*k^2*b_ps*b_ss - (2*k^2 - k_s^2)^2 - (rho_f/rho_s)*k_s^4*(b_ps/b_f) = 0
%
%     where:
%       k    = w / c
%       k_s  = w / vs_s
%       b_ps = sqrt(k^2 - w^2 / vp_s^2)
%       b_ss = sqrt(k^2 - w^2 / vs_s^2)
%       b_f  = sqrt(k^2 - w^2 / vf^2)
%
%     Scattering Cross-Sections (Snieder 1986, Xu et al. 2021 extension):
%       Modal normalization:
%         Inorm = int_{-inf}^0 (|ux_s|^2 + |uz_s|^2) dz + (rho_f/rho_s) int_0^inf (|ux_f|^2 + |uz_f|^2) dz
%
%       In-plane 2D scattering (Scholte -> Scholte):
%         sig_SS = w^3 / (8 * c_S^3 * Inorm^2) * [ Ir4_s * e2_s + (rho_f/rho_s)^2 * Ir4_f * e2_f ]
%
%       Bulk leakage radiation:
%         Scholte -> Fluid Acoustic P_f:
%           sig_SPf = w*kf^3 / (12*pi*c_S) * [(rho_f/rho_s)*Inorm_f / Inorm] * e2_f
%         Scholte -> Solid P:
%           sig_SPs = w*kPs^3 / (12*pi*c_S) * [Inorm_s / Inorm] * e2_s
%         Scholte -> Solid Sv:
%           sig_SSvs = w*kSs^3 / (24*pi*c_S) * [Ir14r2_s / Inorm] * e2_s
%         Scholte -> Solid Sh:
%           sig_SShs = w*kSs^3 / (8*pi*c_S) * [Ir1_s / Inorm] * e2_s
%
%       Total cross section:
%         Sigma_total = n_s*(sig_SS_s + sig_SPs + sig_SSvs + sig_SShs) + n_f*(sig_SS_f + sig_SPf)
%         meanFreePath = 1 / Sigma_total
%
%   Inputs:
%     vp_s, vs_s, rho_s : P-wave, S-wave speed (m/s) and density (kg/m^3) of solid
%     vf, rho_f         : Sound speed (m/s) and density of fluid (vs=0)
%     freq              : Frequency in Hz (default: 1000)
%     het_s, het_f      : Optional heterogeneity structs with fields:
%                         .V (volume, m^3), .drho (delta_rho/rho), .n (density, m^-3)
%
%   Outputs in struct M:
%     .c_Scholte     : Phase velocity (m/s)
%     .deltaFluid    : 1-(c_Scholte/vf)^2 for a fluid-limited mode
%     .cDeficit      : Stable representation of vf-c_Scholte (m)
%     .k_Scholte     : Wavenumber w/c (rad/m)
%     .lambda_S      : Wavelength (m)
%     .b_ps, .b_ss   : Decay rates in solid (1/m)
%     .b_f           : Decay rate in fluid (1/m)
%     .u_x_solid     : Handle for horizontal displacement in solid
%     .u_z_solid     : Handle for vertical displacement in solid
%     .u_x_fluid     : Handle for horizontal displacement in fluid
%     .u_z_fluid     : Handle for vertical displacement in fluid
%     .W_solid       : Energy density in solid (J/m^3)
%     .W_fluid       : Energy density in fluid (J/m^3)
%     .E_solid       : Integrated energy in solid
%     .E_fluid       : Integrated energy in fluid
%     .frac_solid    : Energy fraction in solid (eta_s)
%     .frac_fluid    : Energy fraction in fluid (eta_f)
%     .Inorm         : Modal energy normalization integral
%     .Inorm_solid, .Inorm_fluid : Per-medium parts of Inorm (fluid unweighted)
%     .U_Scholte     : Group velocity (= phase velocity, non-dispersive here)
%     .c_Ps, .c_Ss   : Bulk-to-Scholte coefficients in the solid (m^4)
%     .c_Pf          : Acoustic-to-Scholte coefficient in the fluid (m^4)
%     .sigBar_PsS, .sigBar_SvsS, .sigBar_ShsS, .sigBar_PfS :
%                      angle-averaged, depth-integrated bulk-to-Scholte
%                      cross sections (m^3); see the header note on units
%     .lameSolid, .lameFluid : delta-lambda and delta-mu channels, both
%                      directions, from interfaceLameChannels. Set
%                      het.dlam / het.dmu to activate them; they are zero
%                      by default and then contribute nothing to Sigma.
%     .sig_SS        : In-plane Scholte-to-Scholte cross section (m^2)
%     .sig_SS_solid, .sig_SS_fluid : Per-medium parts of sig_SS (m^2)
%     .sig_SPf       : Scholte-to-fluid-P cross section (m^2)
%     .sig_SPs       : Scholte-to-solid-P cross section (m^2)
%     .sig_SSvs      : Scholte-to-solid-Sv cross section (m^2)
%     .sig_SShs      : Scholte-to-solid-Sh cross section (m^2)
%     .Sigma_solid, .Sigma_fluid : Per-medium scattering rates (1/m)
%     .Sigma_total   : Total scattering rate (1/m)
%     .meanFreePath  : Mean free path ell = 1/Sigma_total (m)
%     .meanFreeTime  : Mean free time tau = ell / c_Scholte (s)
%
%   References:
%     - Scholte, J.G. (1947), Mon. Not. R. Astron. Soc., Geophys. Suppl. 5, 120-126.
%       doi:10.1111/j.1365-246X.1947.tb00347.x
%     - Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Geophys. J. Int. 228, 1213-1236.
%       doi:10.1093/gji/ggab370
%     - Snieder, R. (1986), Geophys. J. Int. 84, 581-605. doi:10.1111/j.1365-246X.1986.tb04372.x

if nargin < 6 || isempty(freq), freq = 1000; end
validateScholteInputs(vp_s, vs_s, rho_s, vf, rho_f, freq);
w = 2 * pi * freq;

% Parse heterogeneity parameters (default: 5% density fluctuations, unit scatterer)
if nargin < 7 || isempty(het_s)
    het_s = struct('V', 1.0, 'drho', 0.05, 'n', 1.0);
elseif isnumeric(het_s)
    het_s = struct('V', sqrt(het_s), 'drho', 1.0, 'n', 1.0);
end
if nargin < 8 || isempty(het_f)
    het_f = struct('V', 1.0, 'drho', 0.05, 'n', 1.0);
elseif isnumeric(het_f)
    het_f = struct('V', sqrt(het_f), 'drho', 1.0, 'n', 1.0);
end
het_s = withLameDefaults(het_s);
het_f = withLameDefaults(het_f);
validateHeterogeneity(het_s, 'het_s');
validateHeterogeneity(het_f, 'het_f');

e2_s = het_s.V^2 * het_s.drho^2;
e2_f = het_f.V^2 * het_f.drho^2;
n_s  = het_s.n;
n_f  = het_f.n;

% A non-zero density fluid always requires a bound Scholte mode below both
% vf and vs_s.  Only the exact rho_f = 0 case is the Rayleigh/vacuum limit.
if rho_f == 0
    c_limit = vs_s;
else
    c_limit = min(vf, vs_s);
end

% When the fluid sound speed is the limiting velocity, solve in
% deltaFluid = 1-(c/vf)^2.  This retains roots arbitrarily close to vf even
% when their velocity deficit is smaller than one double-precision ulp.
fluid_limited = rho_f > 0 && vf <= vs_s;
deltaFluid = NaN;
cDeficit = NaN;

% Keep clear of c = 0 and c = c_limit, where the dimensionless secular
% equation respectively loses precision and can have a square-root
% singularity.  These are floating-point, not physical, cut-offs.
c_min = c_limit * 1e-12;
c_max = c_limit - 8 * eps(c_limit);

% Scholte secular equation function.  It is dimensionless and only defined
% on the physical bound-mode interval.
secular = @(c) computeSecular(c, w, vp_s, vs_s, rho_s, vf, rho_f);

if fluid_limited
    [deltaFluid, c_Scholte, cDeficit] = solveFluidLimitedScholte( ...
        vp_s, vs_s, rho_s, vf, rho_f);
else
    % Use a logarithmic grid in the distance from both endpoints.  A uniform
    % velocity grid misses slow roots in strongly loaded interfaces.
    c_grid = scholteRootGrid(c_min, c_max, c_limit);
    f_grid = arrayfun(secular, c_grid);
    finite_real = isfinite(f_grid) & isreal(f_grid);
    sign_changes = find(finite_real(1:end-1) & finite_real(2:end) & ...
        f_grid(1:end-1) .* f_grid(2:end) <= 0);

    if isempty(sign_changes)
        if rho_f > 0
            error('scholteModel:UnresolvedRoot', ...
                'Could not resolve a bound Scholte root in double precision.');
        end
        error('scholteModel:RootNotBracketed', 'Could not bracket the Rayleigh-limit root.');
    end
    idx = sign_changes(end);
    c_Scholte = fzero(secular, [c_grid(idx), c_grid(idx+1)]);
end

if fluid_limited
    valid_root = isreal(c_Scholte) && c_Scholte > 0 && ...
        isfinite(deltaFluid) && deltaFluid > 0 && deltaFluid < 1;
else
    valid_root = isreal(c_Scholte) && c_Scholte > 0 && c_Scholte < c_limit;
end
if ~valid_root
    error('scholteModel:InvalidRoot', 'The computed phase velocity is not a bound interface-wave root.');
end

ks = w / vs_s;
if fluid_limited
    c_ratio = sqrt(1 - deltaFluid);
    k = w / (vf * c_ratio);
    b_ps = k * sqrt(1 - (vf/vp_s)^2 * (1 - deltaFluid));
    b_ss = k * sqrt(1 - (vf/vs_s)^2 * (1 - deltaFluid));
    b_f = k * sqrt(deltaFluid);
else
    k = w / c_Scholte;
    b_ps = sqrt(k^2 - (w/vp_s)^2);
    b_ss = sqrt(k^2 - (w/vs_s)^2);
end

% Fluid decay factor
if rho_f == 0
    b_f = 0;
    uz_fluid = @(z) zeros(size(z));
    ux_fluid = @(z) zeros(size(z));
    Wf = @(z) zeros(size(z));
    E_fluid = 0.0;
else
    if ~fluid_limited
        b_f = sqrt(k^2 - (w/vf)^2);
    end
    uz_fluid = @(z) exp(-b_f * z);
    ux_fluid = @(z) -1i * (k / b_f) * exp(-b_f * z);

    % Time-averaged acoustic energy: kinetic plus compressional energy.
    % The displacement potential is phi = -exp(-b_f*z)/b_f, hence
    % p = rho_f*w^2*phi.  The sign is immaterial for energy.
    p_fluid = @(z) -rho_f * w^2 / b_f * exp(-b_f * z);
    Kf = @(z) 0.25 * rho_f * w^2 * ...
        (abs(ux_fluid(z)).^2 + abs(uz_fluid(z)).^2);
    Uf = @(z) abs(p_fluid(z)).^2 / (4 * rho_f * vf^2);
    Wf = @(z) Kf(z) + Uf(z);
    E_fluid = rho_f * w^2 * (1.0 + (k/b_f)^2) / (8 * b_f) + ...
               rho_f * w^4 / (8 * vf^2 * b_f^3);
end

% Solid displacement eigenfunctions (normalized to unit vertical displacement at z=0)
denom = 2*k^2 - ks^2;
uz0_raw = b_ps * (1.0 - (2*k^2 / denom));
A = 1.0 / uz0_raw;

ux_solid = @(z) 1i * k * A * (exp(b_ps*z) - (2*b_ps*b_ss/denom) * exp(b_ss*z));
uz_solid = @(z) A * b_ps * (exp(b_ps*z) - (2*k^2/denom) * exp(b_ss*z));

dux_solid_dz = @(z) 1i * k * A * ...
    (b_ps * exp(b_ps*z) - (2*b_ps*b_ss^2/denom) * exp(b_ss*z));
duz_solid_dz = @(z) A * b_ps * ...
    (b_ps * exp(b_ps*z) - (2*k^2*b_ss/denom) * exp(b_ss*z));
Ws = elasticEnergyDensity(ux_solid, uz_solid, dux_solid_dz, duz_solid_dz, ...
    k, rho_s, vp_s, vs_s, w);

% Every eigenfunction is a sum of decaying exponentials, so the semi-infinite
% depth integrals converge and integral() transforms them itself.  A fixed
% 30/min(b) window is both more code and less accurate: for a weakly bound
% mode it spans ~1e8 m and the adaptive quadrature then drifts 2.6e-9 from
% the exact closed form, against 2e-14 here.
E_solid = integral(Ws, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);

E_total = E_solid + E_fluid;
frac_solid = E_solid / E_total;
frac_fluid = E_fluid / E_total;

% ---- Exact Modal Depth Integrals for Scattering Cross-Sections ----
r1_s = @(z) abs(ux_solid(z));
r2_s = @(z) abs(uz_solid(z));
r1_f = @(z) abs(ux_fluid(z));
r2_f = @(z) abs(uz_fluid(z));

% Inorm = int (r1^2 + r2^2) dz
Inorm_s = integral(@(z) r1_s(z).^2 + r2_s(z).^2, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);
if rho_f > 0 && b_f > 0
    Inorm_f = (1.0 + (k/b_f)^2) / (2 * b_f);
    Inorm   = Inorm_s + (rho_f / rho_s) * Inorm_f;
else
    Inorm_f = 0.0;
    Inorm   = Inorm_s;
end

% Quartic integrals Ir4 = int (r1^4 + 2*r2^4) dz
Ir4_s = integral(@(z) r1_s(z).^4 + 2*r2_s(z).^4, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);
if rho_f > 0 && b_f > 0
    Ir4_f = ( (k/b_f)^4 + 2.0 ) / (4 * b_f);
else
    Ir4_f = 0.0;
end

% Quadratic components in solid
Ir1_s    = integral(@(z) r1_s(z).^2, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);
Ir2_s    = integral(@(z) r2_s(z).^2, -Inf, 0, 'RelTol', 1e-10, 'AbsTol', 0);
Ir14r2_s = Ir1_s + 4 * Ir2_s;

% ---- Cross-Section Calculations ----
kPs = w / vp_s;
kSs = w / vs_s;
if rho_f > 0 && vf > 0
    kPf = w / vf;
else
    kPf = 0.0;
end

% 1. In-plane 2D scattering (Scholte -> Scholte)
sig_SS_s = (w^3 / (8 * c_Scholte^3 * Inorm^2)) * Ir4_s * e2_s;
if rho_f > 0 && b_f > 0
    sig_SS_f = (w^3 / (8 * c_Scholte^3 * Inorm^2)) * ((rho_f/rho_s)^2 * Ir4_f) * e2_f;
else
    sig_SS_f = 0.0;
end
sig_SS = sig_SS_s + sig_SS_f;

% 2. Bulk radiation into fluid acoustic wave (Scholte -> P_fluid)
if rho_f > 0 && b_f > 0
    sig_SPf = (w * kPf^3 / (12 * pi * c_Scholte)) * ((rho_f/rho_s) * Inorm_f / Inorm) * e2_f;
else
    sig_SPf = 0.0;
end

% 3. Bulk radiation into solid (Scholte -> P_solid, Sv_solid, Sh_solid)
% The P kernel is (r1^2 + r2^2) (Xu et al., 2021, B.6), so its depth integral
% is the per-medium part of Inorm.  In a single halfspace that ratio is 1,
% which is why B.17 shows none; here the two media share one normalization,
% so Inorm_solid/Inorm and (rho_f/rho_s)*Inorm_fluid/Inorm must sum to 1.
sig_SPs  = (w * kPs^3 / (12 * pi * c_Scholte)) * (Inorm_s / Inorm) * e2_s;
sig_SSvs = (w * kSs^3 / (24 * pi * c_Scholte)) * (Ir14r2_s / Inorm) * e2_s;
sig_SShs = (w * kSs^3 / ( 8 * pi * c_Scholte)) * (Ir1_s / Inorm) * e2_s;

% ---- Inverse channels: bulk -> Scholte -------------------------------
% The reverse of the three solid channels and the fluid one, i.e. the
% analogue of RayleighModel's cP/cS and Xu et al. B.3-B.5. A bulk wave of
% mode b arriving at incidence psi from the interface normal, scattering off
% a density heterogeneity at depth z, converts into Scholte with
%
%   sigma_{P ->S}(z,psi) = c_P*[ r1(z)^2 sin^2(psi) + 2 r2(z)^2 cos^2(psi) ]
%   sigma_{Sv->S}(z,psi) = c_S*[ r1(z)^2 cos^2(psi) + 2 r2(z)^2 sin^2(psi) ]
%   sigma_{Sh->S}(z)     = c_S*  r1(z)^2
%
% and the fluid acoustic wave uses the P kernel, because e_A = n_hat is
% longitudinal just like e_P. Units: c_* are m^4, the kernels m^-2, so each
% sigma is m^2 AT A GIVEN DEPTH -- unlike the guided->bulk channels above,
% which are already integrated over the scatterer depth. The two directions
% are not interchangeable and are not equal; see the density-of-states ratio
% asserted in testInterfaceWavesAndWedge.m.
%
% The (rho_j/rho_s) weight is the same one the forward channels carry, and it
% has to be: it cancels in detailed balance, leaving the medium-independent
% ratio pi*v^3/(w*c*U) that RayleighModel already satisfies to 2e-16.
U_S = c_Scholte;                        % no dispersion between two half-spaces
c_Ps = (w^3 / (8 * vp_s * c_Scholte * U_S * Inorm)) * e2_s;
c_Ss = (w^3 / (8 * vs_s * c_Scholte * U_S * Inorm)) * e2_s;
if rho_f > 0 && b_f > 0
    c_Pf = (w^3 / (8 * vf * c_Scholte * U_S * Inorm)) * (rho_f/rho_s) * e2_f;
else
    c_Pf = 0.0;
end

% Angle-averaged and depth-integrated versions (the B.6-B.8 analogue), which
% are what detailed balance compares against.
sigBar_PsS  = (2/3) * c_Ps * Inorm_s;
sigBar_SvsS = (1/3) * c_Ss * (Ir1_s + 4*Ir2_s);
sigBar_ShsS =         c_Ss * Ir1_s;
sigBar_PfS  = (2/3) * c_Pf * Inorm_f;

% ---- delta-lambda and delta-mu channels -------------------------------
% Same two directions again, now for heterogeneity in the Lame parameters.
% Zero unless het.dlam / het.dmu were set, so nothing an existing caller
% sees can change. The fluid has mu = 0, hence delta-lambda only, into the
% acoustic wave.
Ks_lame = interfaceLameKernels(ux_solid, uz_solid, dux_solid_dz, duz_solid_dz, k);
Cs_lame = interfaceLameChannels(Ks_lame, [-Inf 0], w, vp_s, vs_s, rho_s, rho_s, ...
    c_Scholte, U_S, Inorm, het_s.V, het_s.dlam, het_s.dmu);
if rho_f > 0 && b_f > 0
    dux_fluid_dz = @(z) 1i * k * exp(-b_f * z);
    duz_fluid_dz = @(z) -b_f * exp(-b_f * z);
    Kf_lame = interfaceLameKernels(ux_fluid, uz_fluid, dux_fluid_dz, duz_fluid_dz, k);
    Cf_lame = interfaceLameChannels(Kf_lame, [0 Inf], w, vf, 0, rho_f, rho_s, ...
        c_Scholte, U_S, Inorm, het_f.V, het_f.dlam, 0);
else
    Cf_lame = interfaceLameChannels(Ks_lame, [-Inf 0], w, vp_s, vs_s, rho_s, rho_s, ...
        c_Scholte, U_S, Inorm, 0, 0, 0);
end

% Total cross section Sigma (1/m)
Sigma_solid_channels = n_s * (sig_SS_s + sig_SPs + sig_SSvs + sig_SShs ...
    + Cs_lame.sig_lam_P + Cs_lame.sig_mu_P + Cs_lame.sig_mu_Sv + Cs_lame.sig_mu_Sh);
Sigma_fluid_channels = n_f * (sig_SS_f + sig_SPf + Cf_lame.sig_lam_P);
Sigma_total = Sigma_solid_channels + Sigma_fluid_channels;

meanFreePath = 1.0 / Sigma_total;
meanFreeTime = meanFreePath / c_Scholte;

% Effective scattering proxy:
eff_sigma = @(sig_s, sig_f) frac_solid * sig_s + frac_fluid * sig_f;

% Populate struct
M.c_Scholte     = c_Scholte;
M.deltaFluid    = deltaFluid;
M.cDeficit      = cDeficit;
M.k_Scholte     = k;
M.lambda_S      = 2*pi / k;
M.b_ps          = b_ps;
M.b_ss          = b_ss;
M.b_f           = b_f;
M.u_x_solid     = ux_solid;
M.u_z_solid     = uz_solid;
M.u_x_fluid     = ux_fluid;
M.u_z_fluid     = uz_fluid;
M.W_solid       = Ws;
M.W_fluid       = Wf;
M.E_solid       = E_solid;
M.E_fluid       = E_fluid;
M.frac_solid    = frac_solid;
M.frac_fluid    = frac_fluid;
M.Inorm         = Inorm;
M.Inorm_solid   = Inorm_s;
M.Inorm_fluid   = Inorm_f;
M.U_Scholte     = U_S;
M.c_Ps          = c_Ps;
M.c_Ss          = c_Ss;
M.c_Pf          = c_Pf;
M.sigBar_PsS    = sigBar_PsS;
M.sigBar_SvsS   = sigBar_SvsS;
M.sigBar_ShsS   = sigBar_ShsS;
M.sigBar_PfS    = sigBar_PfS;
M.lameSolid     = Cs_lame;
M.lameFluid     = Cf_lame;
M.sig_SS        = sig_SS;
M.sig_SS_solid  = sig_SS_s;
M.sig_SS_fluid  = sig_SS_f;
M.sig_SPf       = sig_SPf;
M.sig_SPs       = sig_SPs;
M.sig_SSvs      = sig_SSvs;
M.sig_SShs      = sig_SShs;
M.Sigma_solid   = Sigma_solid_channels;
M.Sigma_fluid   = Sigma_fluid_channels;
M.Sigma_total   = Sigma_total;
M.meanFreePath  = meanFreePath;
M.meanFreeTime  = meanFreeTime;
M.eff_sigma     = eff_sigma;

end

function val = computeSecular(c, ~, vp_s, vs_s, rho_s, vf, rho_f)
    if ~isreal(c) || c <= 0
        val = NaN;
        return;
    end

    if rho_f == 0
        c_limit = vs_s;
    else
        c_limit = min(vf, vs_s);
    end
    if c >= c_limit
        val = NaN;
        return;
    end

    y = (c / vs_s)^2;
    gamma2 = (vs_s / vp_s)^2;
    ap = sqrt(1 - gamma2 * y);
    as = sqrt(1 - y);

    % Evaluate the Rayleigh part without subtracting two numbers that both
    % approach four as c -> 0.  Dividing by y removes the trivial c = 0
    % zero of the original dimensional determinant without changing any
    % finite-speed root.
    A = ap * as;
    B2 = (1 - y / 2)^2;
    rayleigh_numerator = (1 - gamma2) * y + (gamma2 - 1.5) * y^2 + ...
                         0.5 * y^3 - y^4 / 16;
    rayleigh_part = 4 * rayleigh_numerator / (A + B2);

    if rho_f == 0
        fluid_part = 0;
    else
        af = sqrt(1 - (c / vf)^2);
        fluid_part = (rho_f / rho_s) * y^2 * (ap / af);
    end

    val = (rayleigh_part - fluid_part) / y;
end

function [delta, c, c_deficit] = solveFluidLimitedScholte(vp_s, vs_s, rho_s, vf, rho_f)
    % eta = log(delta) gives uniform relative resolution for a root that may
    % be exponentially close to the fluid sound speed.
    eta_min = log(realmin);
    eta_max = log(1 - 1e-12);
    secular_eta = @(eta) computeSecularFromDelta(exp(eta), vp_s, vs_s, rho_s, vf, rho_f);
    eta_grid = linspace(eta_min, eta_max, 1201);
    values = arrayfun(secular_eta, eta_grid);
    finite_real = isfinite(values) & isreal(values);
    sign_changes = find(finite_real(1:end-1) & finite_real(2:end) & ...
        values(1:end-1) .* values(2:end) <= 0);
    if isempty(sign_changes)
        error('scholteModel:UnresolvedRoot', ...
            'Could not bracket the fluid-limited Scholte root in log(delta).');
    end

    % The physical fluid-loaded Scholte root is the one closest to vf, i.e.
    % the smallest positive delta.
    idx = sign_changes(1);
    eta = fzero(secular_eta, [eta_grid(idx), eta_grid(idx + 1)]);
    delta = exp(eta);
    c_ratio = sqrt(1 - delta);
    c = vf * c_ratio;
    % This rationalized expression remains meaningful when c rounds to vf.
    c_deficit = vf * delta / (1 + c_ratio);
end

function val = computeSecularFromDelta(delta, vp_s, vs_s, rho_s, vf, rho_f)
    if ~isreal(delta) || delta <= 0 || delta >= 1
        val = NaN;
        return;
    end

    c_ratio_sq = 1 - delta;
    y = (vf / vs_s)^2 * c_ratio_sq;
    gamma2 = (vs_s / vp_s)^2;
    ap = sqrt(1 - gamma2 * y);
    as = sqrt(1 - y);
    A = ap * as;
    B2 = (1 - y / 2)^2;
    rayleigh_numerator = (1 - gamma2) * y + (gamma2 - 1.5) * y^2 + ...
                         0.5 * y^3 - y^4 / 16;
    rayleigh_part = 4 * rayleigh_numerator / (A + B2);
    fluid_part = (rho_f / rho_s) * y^2 * (ap / sqrt(delta));
    val = (rayleigh_part - fluid_part) / y;
end

function c_grid = scholteRootGrid(c_min, c_max, c_limit)
    % Cover roots near zero with log(c/c_limit), and roots near c_limit with
    % log((c_limit-c)/c_limit).  Endpoints stay strictly inside the domain.
    q_min = c_min / c_limit;
    q_max = c_max / c_limit;
    d_min = 1 - q_max;

    q_low  = logspace(log10(q_min), -1, 350);
    q_mid  = linspace(0.1, 0.9, 401);
    d_high = logspace(log10(d_min), -1, 350);
    q_high = 1 - fliplr(d_high);

    q = unique([q_low, q_mid, q_high]);
    q = q(q >= q_min & q <= q_max);
    c_grid = c_limit * q;
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

function validateScholteInputs(vp_s, vs_s, rho_s, vf, rho_f, freq)
    validatePositiveScalar(vp_s, 'vp_s');
    validatePositiveScalar(vs_s, 'vs_s');
    validatePositiveScalar(rho_s, 'rho_s');
    validatePositiveScalar(vf, 'vf');
    validateNonnegativeScalar(rho_f, 'rho_f');
    validatePositiveScalar(freq, 'freq');
    if vp_s <= vs_s
        error('scholteModel:InvalidInput', 'vp_s must be greater than vs_s for an isotropic elastic solid.');
    end
end

function het = withLameDefaults(het)
    % Lame perturbations are optional and default to zero, so every existing
    % caller keeps exactly the cross sections it had.
    if ~isfield(het, 'dlam') || isempty(het.dlam), het.dlam = 0; end
    if ~isfield(het, 'dmu')  || isempty(het.dmu),  het.dmu  = 0; end
end

function validateHeterogeneity(het, name)
    if ~isstruct(het) || ~isscalar(het) || ~all(isfield(het, {'V', 'drho', 'n'}))
        error('scholteModel:InvalidHeterogeneity', ...
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
        error('scholteModel:InvalidInput', '%s must be a positive finite real scalar.', name);
    end
end

function validateNonnegativeScalar(value, name)
    validateFiniteRealScalar(value, name);
    if value < 0
        error('scholteModel:InvalidInput', '%s must be a non-negative finite real scalar.', name);
    end
end

function validateFiniteRealScalar(value, name)
    if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value)
        error('scholteModel:InvalidInput', '%s must be a finite real scalar.', name);
    end
end
