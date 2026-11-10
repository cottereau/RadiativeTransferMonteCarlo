# MATLAB Architecture — 3D Extended Lamb's Problem, DWN + Energy Partition

Planned signatures, I/O structures and docstrings (§1–§7). **Implemented** in `matlab/+lamb` — see §8 for the deviations.
Equation numbers `(eq:…)` refer to labels in `LambFK_Formulation.tex`.

---

## 0. Fixed conventions (all modules)

| Item | Convention |
|---|---|
| Time dependence | `exp(-1i*omega*t)`; complex frequency `omegaC = omega + 1i*eps` (eps > 0) |
| Geometry | interface z = 0, z positive **down**; solid (medium 1) in z > 0 |
| Source | at (0,0,zs); `zs > 0` → in solid; `zs < 0` → in upper medium (fluid, Case 3) |
| Wavenumbers | polar (k, phi); local frame l̂, t̂, ẑ; arrays are `Nk × Nphi` |
| Vertical wavenumber | `nu = -1i*sqrt(km.^2 - k.^2)` → Re ν ≥ 0 (valid for real and complex ω) |
| Amplitudes | P, SV = potentials (φ, ψ); SH = displacement u_T |
| Units | SI inside the classes; solvers internally nondimensionalise by ρ₁, β₁, 1/k_β₁ |
| MATLAB | R2022a+ (`pagemldivide`, `pagemtimes`); no toolboxes required |

---

## 1. Folder layout

```
RadiativeTransferMonteCarlo/lamb/
├── +lamb/
│   ├── Medium.m                 % classdef (value): elastic solid or acoustic fluid
│   ├── Interface.m              % classdef (value): lower solid + upper medium / free surface
│   ├── PointSource.m            % classdef (value): force (any direction) or moment tensor
│   ├── incidentAmplitudes.m     % (eq:PSH) → P0,S0,C0 (toward interface) + Pd0,Sd0,Cd0 (away)
│   ├── solveInterface.m         % M x = S, vectorised over (k,phi) pages
│   ├── findPoles.m              % Rayleigh / Scholte / Stoneley roots on real k axis
│   ├── poleResidues.m           % (eq:res) residues via null vectors
│   ├── energyPartition.m        % *** main energy module (driver) ***
│   ├── +energy/
│   │   ├── fluxDensity.m        % E_m(k) (eq:Em) and E_SH
│   │   ├── bodyWaves.m          % (eq:Wbody3)
│   │   ├── interfaceWaves.m     % (eq:WIW3) flux + source-side value
│   │   ├── sourcePower.m        % (eq:Wsrc3) independent balance check
│   │   ├── rtCoefficients.m     % plane-wave R/T energy coefficients vs slowness
│   │   └── broadband.m          % wavelet-weighted partition over ω
│   ├── +quad/
│   │   ├── windowNodes.m        % branch-point-split Gauss nodes (eq:subst)
│   │   └── azimuthNodes.m       % exact uniform trapezoid in phi
│   └── +wavefield/
│       ├── dwn.m                % DWN spectra U(ω_c) at receivers (azimuthal FFT + Bessel sum)
│       └── synthetics.m         % time series: complex-ω fft synthesis, exp(eps t) correction
├── tests/
│   ├── testWholeSpace.m
│   ├── testMillerPursey.m
│   ├── testEnergyBalance.m
│   ├── testRTUnitarity.m
│   ├── testFreeSurfaceLimit.m
│   └── testSourceSideSymmetry.m
└── examples/
    └── runPartitionVsDepth.m
```

Design choices: value classes (no handles, no hidden state), and all heavy functions
take `(iface, src, …)` explicitly. Everything is vectorised over `(k, phi)` using `Nk × Nphi`
arrays and `4×4×(Nk·Nphi)` page matrices. There is no problem/model god-object; a plain
struct can bundle `iface` and `src` in scripts.

---

## 2. Classes

### 2.1 `lamb.Medium`

```matlab
classdef Medium
    %MEDIUM  Homogeneous isotropic elastic solid or inviscid acoustic fluid.
    %
    %   m = lamb.Medium(rho, alpha, beta)   % solid (beta > 0)
    %   m = lamb.Medium(rho, alpha, 0)      % fluid (beta = 0)
    %
    %   Properties
    %     rho    density [kg/m^3]
    %     alpha  P (sound) speed [m/s]
    %     beta   S speed [m/s]; 0 for a fluid
    %   Dependent
    %     mu, lambda   Lamé parameters (mu = 0 for a fluid)
    %     isFluid      logical
    %
    %   Validation: rho > 0, alpha > 0, 0 <= beta < alpha/sqrt(2) (solid).
    properties (SetAccess = immutable)
        rho   (1,1) double {mustBePositive}
        alpha (1,1) double {mustBePositive}
        beta  (1,1) double {mustBeNonnegative}
    end
    properties (Dependent)
        mu; lambda; isFluid
    end
    methods
        function obj = Medium(rho, alpha, beta)
        end
        function nu = verticalWavenumber(obj, k, omega, wave)
            %VERTICALWAVENUMBER  nu_m = sqrt(k^2 - (omega/c_m)^2), Re >= 0 (eq:branch).
            %   nu = m.verticalWavenumber(k, omega, 'P' | 'S')
            %   k      Nk×Nphi (or any size) real wavenumbers
            %   omega  scalar, real or complex (omega + 1i*eps)
            %   nu     same size as k, complex
        end
        function km = wavenumber(obj, omega, wave)
            %WAVENUMBER  k_m = omega / c_m for wave = 'P' | 'S'.
        end
    end
end
```

### 2.2 `lamb.Interface`

```matlab
classdef Interface
    %INTERFACE  Welded plane interface z = 0 between a lower solid and an upper medium.
    %
    %   I = lamb.Interface(lowerSolid, upperMedium)
    %   I = lamb.Interface(lowerSolid, [])            % traction-free surface (Lamb)
    %
    %   Properties
    %     lower   lamb.Medium, must be solid (occupies z > 0)
    %     upper   lamb.Medium or [] (occupies z < 0)
    %   Dependent
    %     kind    'solid-solid' | 'solid-fluid' | 'free'
    %
    %   The P-SV matrix depends only on the interface, never on the source:
    %     'solid-solid' -> M1 (4x4, eq:case1),   SH 2x2
    %     'solid-fluid' -> M2 (3x3, eq:case2),   SH: C1 = C0
    %     'free'        -> M_FS (2x2),           SH: C1 = C0
    properties (SetAccess = immutable)
        lower lamb.Medium
        upper            % lamb.Medium or []
    end
    properties (Dependent)
        kind
    end
    methods
        function obj = Interface(lower, upper)
        end
        function [M, Msh] = matrices(obj, k, omega)
            %MATRICES  Boundary-condition matrices as pages.
            %   [M, Msh] = I.matrices(k, omega)
            %   k      N×1 (or Nk×Nphi, reshaped internally)
            %   M      n×n×N complex, n = 4 | 3 | 2 (by kind); stress rows scaled by 1/mu1
            %   Msh    2×2×N (solid-solid) or [] (SH decoupled trivially)
        end
        function [Delta, dDelta] = dispersion(obj, k, omega)
            %DISPERSION  Delta(k) = det M(k), real-valued (to a constant phase) for k > max branch pt.
            %   Also returns dDelta/dk (analytic). Scaled so that Delta = O(1).
        end
        function kb = branchPoints(obj, omega)
            %BRANCHPOINTS  Sorted real branch points {k_alpha1, k_beta1, k_alpha2, k_beta2, k_f}.
        end
    end
end
```

### 2.3 `lamb.PointSource`

```matlab
classdef PointSource
    %POINTSOURCE  Point force of arbitrary direction or arbitrary moment tensor at (0,0,zs).
    %
    %   s = lamb.PointSource('force',  F,  zs)    % F  3×1 [N]  (x, y, z-down)
    %   s = lamb.PointSource('moment', M,  zs)    % M  3×3 [N·m], any matrix (symmetric = seismic)
    %   s = lamb.PointSource(..., 'Wavelet', @(omega) ...)   % source spectrum s(omega), default 1
    %
    %   Helpers (static):
    %     lamb.PointSource.explosion(M0, zs)                 % M = M0*eye(3)
    %     lamb.PointSource.doubleCouple(M0, strike, dip, rake, zs)
    %
    %   zs > 0: source in the lower solid. zs < 0: source in the upper medium (fluid only).
    properties (SetAccess = immutable)
        kind  (1,:) char {mustBeMember(kind, {'force','moment'})} = 'force'
        F     (3,1) double = [0;0;1]
        M     (3,3) double = zeros(3)
        zs    (1,1) double
        wavelet function_handle = @(omega) ones(size(omega))
    end
    methods
        function obj = PointSource(kind, FM, zs, varargin)
        end
        function [Pi, Sig, H] = directAmplitudes(obj, med, k, phi, omega, epsSign)
            %DIRECTAMPLITUDES  Whole-space up/down-going amplitudes (eq:srcF, eq:srcM, eq:PSH).
            %   med      lamb.Medium containing the source
            %   k, phi   Nk×Nphi
            %   epsSign  +1 (field below source, z > zs) | -1 (above)
            %   Pi       P potential amplitude   of e^{-nu_a|z-zs|}
            %   Sig      SV potential amplitude  of e^{-nu_b|z-zs|}   (0 in a fluid)
            %   H        SH displacement ampl.   of e^{-nu_b|z-zs|}   (0 in a fluid)
        end
        function S = functional(obj, v, kappa)
            %FUNCTIONAL  Source functional S[u] used in the source-work check.
            %   force:  S = F' * u(xs);   moment: S = sum_waves v.' * conj(M) * kappa
            %   v, kappa  3×N polarisation and wave vectors of scattered waves at z = zs.
        end
    end
end
```

---

## 3. Solver functions

```matlab
function amp = incidentAmplitudes(src, iface, k, phi, omega)
%INCIDENTAMPLITUDES  Source field decomposed at the interface (eq:PSH).
%   amp = lamb.incidentAmplitudes(src, iface, k, phi, omega)
%
%   Inputs   k, phi   Nk×Nphi; omega real or complex
%   Output   amp struct (all Nk×Nphi complex):
%     source in solid (zs = h > 0):
%       .P0 .S0 .C0      toward interface, ∝ e^{+nu z}  (enter S of M x = S)
%       .Pd0 .Sd0 .Cd0   away from interface, ∝ e^{-nu z} (for body flux, |k| < k_m only)
%     source in fluid (zs = -h < 0):
%       .Pf0             toward interface, ∝ e^{-nu_f z}
%       .Uf0             away (upward),    ∝ e^{+nu_f z}
%   Note: Pd0/Sd0/Cd0/Uf0 contain e^{+nu h}; evaluate only inside propagating windows.
```

```matlab
function x = solveInterface(iface, amp, k, omega)
%SOLVEINTERFACE  Solve M x = S (P-SV) and the SH system for all (k,phi) pages.
%   x = lamb.solveInterface(iface, amp, k, omega)
%
%   Output x struct (Nk×Nphi complex; absent fields = []):
%     .A1 .B1   reflected (source in solid) or transmitted (source in fluid) P, SV in solid
%     .A2 .B2   transmitted P, SV into solid 2          (solid-solid)
%     .Af       fluid P (transmitted or reflected)       (solid-fluid)
%     .C1 .C2   SH in solid 1 / solid 2
%     .cond     condition-number estimate per page (diagnostic)
%   Uses pagemldivide on 1/mu1-scaled matrices; for 2×2/3×3 closed-form Cramer is allowed.
```

```matlab
function poles = findPoles(iface, omega, opts)
%FINDPOLES  Real interface-wave poles k_p > max branch point (physical sheet).
%   poles = lamb.findPoles(iface, omega)
%   poles = lamb.findPoles(iface, omega, 'kMaxFactor', 5, 'NScan', 2000)
%
%   Scans Delta(k) for sign changes on (max(branchPoints), kMaxFactor·max(branchPoints)],
%   refines with fzero. Returns struct array (possibly empty, e.g. no Stoneley):
%     .k      real pole wavenumber
%     .c      phase velocity omega/k
%     .type   'Rayleigh' | 'Scholte' | 'Stoneley'
%     .dDelta dDelta/dk at the pole
%   Leaky (complex) poles are NOT returned: their energy is body-wave energy (sec:leaky).
arguments
    iface lamb.Interface
    omega (1,1) double {mustBeReal, mustBePositive}
    opts.kMaxFactor (1,1) double = 5
    opts.NScan (1,1) double = 2000
end
```

```matlab
function res = poleResidues(iface, src, pole, phi, omega)
%POLERESIDUES  Residues of the scattered amplitudes at a real pole (eq:res).
%   res = lamb.poleResidues(iface, src, pole, phi, omega)
%
%   Res x = r (l.' S(k_p)) / (l.' M'(k_p) r), with r, l the right/left null vectors of M(k_p)
%   (from svd), M' analytic or centred difference.
%   Output res struct, fields as in solveInterface (.A1 .B1 .A2 .B2 .Af), each 1×Nphi.
%   SH fields are absent (no SH pole for two half-spaces).
```

---

## 4. Energy Partition module (core deliverable)

### 4.1 Driver

```matlab
function EP = energyPartition(iface, src, omega, opts)
%ENERGYPARTITION  Radiated-power partition between body waves and interface waves.
%   EP = lamb.energyPartition(iface, src, omega)
%   EP = lamb.energyPartition(iface, src, omega, 'NGauss', 64, 'NPhi', 8, ...)
%
%   Evaluated at REAL omega (lossless). Poles are removed analytically by residues;
%   body-wave windows |k| < k_m contain no poles. (The complex-frequency trick is for
%   wavefield synthesis only; with eps > 0 energy is dissipated and the IW/body split blurs.)
%
%   Options
%     NGauss        Gauss-Legendre nodes per branch-point segment (default 64)
%     NPhi          uniform azimuth nodes, >= 5 exact for force & moment (default 8)
%     CheckBalance  compute source work and relative error (default true)
%     RT            also return plane-wave R/T energy coefficients (default false)
%     Pattern       also return dW/dtheta per body wave (default false)
%
%   Output EP struct
%     .omega
%     .W.body       struct of scalars [W]:  .P1 .SV1 .SH1 .P2 .SV2 .SH2 .Pf  (absent -> 0)
%     .W.interface  struct array per pole:  .type .k .c .W (flux, eq:WIW3) .Wsrc (eq:Wsrc3 pole term)
%     .W.bodyTotal  sum of .W.body
%     .W.ifTotal    sum of .W.interface.W
%     .W.total      bodyTotal + ifTotal
%     .pct          same layout as .W.body plus .interface (vector) and .bodyTotal, .ifTotal  [%]
%     .check.Wsrc          total source work (eq:Wsrc3)
%     .check.relErr        |Wsrc - W.total| / W.total   (expect < 1e-6)
%     .check.ifRelErr      |W - Wsrc| per pole          (flux vs source-side residue)
%     .rt           (opt) output of lamb.energy.rtCoefficients for each incident type
%     .pattern      (opt) .theta, .P1, .SV1, ... dW/dtheta
arguments
    iface lamb.Interface
    src   lamb.PointSource
    omega (1,1) double {mustBeReal, mustBePositive}
    opts.NGauss (1,1) double = 64
    opts.NPhi   (1,1) double = 8
    opts.CheckBalance (1,1) logical = true
    opts.RT      (1,1) logical = false
    opts.Pattern (1,1) logical = false
end
```

### 4.2 Sub-functions (`+lamb/+energy`)

```matlab
function E = fluxDensity(med, k, omega, wave)
%FLUXDENSITY  Vertical power density per unit |amplitude|^2 of one plane wave.
%   P, SV (potential amplitude):  E = 0.5*rho*omega^3*real(1i*nu_m)      (eq:Em)
%   SH  (displacement amplitude): E = 0.5*omega*mu*real(1i*nu_beta)
%   Zero outside |k| < k_m automatically.  wave = 'P' | 'SV' | 'SH'.
```

```matlab
function Wb = bodyWaves(iface, src, omega, opts)
%BODYWAVES  Body-wave powers from far-field vertical flux via Parseval (eq:Wbody3).
%   Wb = lamb.energy.bodyWaves(iface, src, omega, opts)
%   For each wave m: integrate E_m |D_m|^2 k over k in [0,k_m] (split at ALL branch points,
%   lamb.quad.windowNodes) and phi (lamb.quad.azimuthNodes).
%   D_m = direct (away) + scattered, coherently (table in sec:energy + SH rows).
%   Output struct .P1 .SV1 .SH1 .P2 .SV2 .SH2 .Pf [W]; optional .pattern (dW/dtheta).
```

```matlab
function Wi = interfaceWaves(iface, src, poles, omega, opts)
%INTERFACEWAVES  Trapped power in each real pole (eq:WIW3), two independent ways.
%   Wi = lamb.energy.interfaceWaves(iface, src, poles, omega, opts)
%   (a) flux:        W = k_p/(2 pi) * int_0^{2pi} W2D+(phi) dphi, with the depth integral of
%                    (omega/2) Im(s_LL u_L* + s_Lz u_z*) done analytically (1/(nu_m + nu_n) terms)
%   (b) source side: W = omega k_p/(8 pi) * int Re S[Res u_sc] dphi
%   Output struct array: .type .k .c .W (a) .Wsrc (b) .eigenfunction (optional, z-profile)
```

```matlab
function Ws = sourcePower(iface, src, poles, omega, opts)
%SOURCEPOWER  Work done by the source, independent of flux bookkeeping (eq:Wsrc3).
%   Ws = W0 (whole-space, Parseval of direct field; analytic for force)
%      + (omega/8 pi^2) PV int int Im S[u_sc] k dk dphi     (regular part, |k| < k_max window)
%      + sum_poles omega k_p/(8 pi) int Re S[Res u_sc] dphi
%   Output struct .total .W0 .regular .pole (vector per pole)
```

```matlab
function rt = rtCoefficients(iface, incident, p, omega)
%RTCOEFFICIENTS  Plane-wave energy reflection/transmission coefficients vs horizontal slowness.
%   rt = lamb.energy.rtCoefficients(iface, 'P1'|'SV1'|'SH1'|'Pf', p, omega)
%   p      vector of horizontal slowness [s/m] (k = omega*p)
%   rt     struct: .p, .(outgoingType) energy fractions (e.g. .P1 .SV1 .P2 .SV2 .Pf .SH1 .SH2),
%          .sum (== 1 for propagating incidence; unitarity diagnostic)
%   e_m = E_m |x_m|^2 / (E_inc |a_inc|^2)  with E from fluxDensity.
```

```matlab
function EPb = broadband(iface, src, omegaVec, opts)
%BROADBAND  Wavelet-weighted energy partition: W_m = int |s(omega)|^2 W_m(omega) domega.
%   EPb = lamb.energy.broadband(iface, src, omegaVec, opts)
%   Calls energyPartition per omega (parfor-ready). Output: .perOmega (struct array of EP),
%   .W (integrated, same layout as EP.W), .pct.
```

### 4.3 Quadrature (`+lamb/+quad`)

```matlab
function [k, wk] = windowNodes(breaks, n)
%WINDOWNODES  Gauss nodes on [0, k_m] split at interior branch points (eq:subst).
%   breaks  sorted [0, k_n1, ..., k_m];  n nodes per segment
%   On [ka,kb]: k^2 = ka^2 + (kb^2-ka^2) sin^2 t  (k = kb sin t if ka = 0), t Gauss-Legendre.
%   Removes the 1/sqrt singularities at both segment ends (Weyl factor, head waves).
```

```matlab
function [phi, wphi] = azimuthNodes(N)
%AZIMUTHNODES  Uniform nodes phi_j = 2 pi j/N, weights 2 pi/N (exact for degree < N).
```

---

## 5. Wavefield module (DWN, complex frequency)

```matlab
function U = dwn(iface, src, rcv, omegaC, opts)
%DWN  Displacement spectra at receivers by the discrete wavenumber method.
%   U = lamb.wavefield.dwn(iface, src, rcv, omegaC, 'L', L, 'Nk', Nk, 'NPhi', 8)
%
%   rcv     struct .x .y .z (Nr×1) [m]
%   omegaC  complex frequencies omega + 1i*eps (Nw×1)
%   L       source periodicity [m]: L > r_max + c_max*T;  k_j = j*2*pi/L, j = 1..Nk
%   Method  (k,phi) spectra -> FFT over phi (|n| <= 3) -> Bessel sum:
%           u = (dk/2pi) sum_n i^n e^{i n theta} sum_j u_n(k_j,z) J_n(k_j r) k_j
%   Output  U  Nr×3×Nw complex (ux, uy, uz)
```

```matlab
function [t, u] = synthetics(iface, src, rcv, topts)
%SYNTHETICS  Time-domain seismograms with the complex-frequency trick.
%   [t, u] = lamb.wavefield.synthetics(iface, src, rcv, 'T', T, 'dt', dt, 'Eps', pi/T)
%   Builds omegaC = omega + 1i*eps, calls dwn, applies src.wavelet, synthesises with the
%   e^{-i omega t} kernel (MATLAB fft, NOT ifft), multiplies by exp(eps*t).
%   Output t (Nt×1), u (Nt×3×Nr) real.
```

---

## 6. Tests (`matlab.unittest`, function-based)

| Test | Setup | Assertion |
|---|---|---|
| `testWholeSpace` | force F=(0.3,−0.5,0.8), ρ=β=ω=1, α=√3 | W0 = ω²‖F‖²(α⁻³+2β⁻³)/(24πρ) (rel 1e-6); S/P = 2(α/β)³; DC S/P = 1.5(α/β)⁵ |
| `testMillerPursey` | free surface, vertical force, h=0, Poisson | P 6.88 %, S 25.75 %, R 67.36 % (±0.02) |
| `testEnergyBalance` | free surface, M=[1 .4 −.7; .4 −.3 .5; −.7 .5 .6], k_βh=0.5 | relErr < 1e-6; P/SV/SH/R = 2.69/8.48/62.02/26.81 % |
| `testRTUnitarity` | all kinds, all incident types, propagating p | `rt.sum == 1` (1e-10) |
| `testFreeSurfaceLimit` | solid–fluid, Z_f = ρ_fα_f → 0 | Rayleigh pole/partition → free-surface values; fluid share → 0 |
| `testSourceSideSymmetry` | any case | `interfaceWaves` (a) vs (b) agree (1e-8) |

---

## 7. Implementation order (Phase 2)

1. `Medium`, `Interface.matrices`/`dispersion`, `findPoles` (+ Rayleigh c_R = 0.9194β).
2. `PointSource.directAmplitudes`, `incidentAmplitudes` → `testWholeSpace`.
3. `solveInterface` (free surface) → `energyPartition` → `testMillerPursey`, `testEnergyBalance`.
4. Solid–fluid (Cases 2, 3), then solid–solid (Case 1); `rtCoefficients`.
5. `wavefield.dwn` / `synthetics` last (energy module does not depend on them).

---

## 8. Implementation status (2026-09-24) — deviations from the plan above

All of §2–§5 is implemented in `matlab/+lamb` and passes `tests/testLamb.m` (7/7) and
`tests/test_rayleigh_energy_limit.m`. Deviations:

| Planned | Implemented | Why |
|---|---|---|
| `solveInterface.m` | `SolidSolidSolver` (Case 1) and `SolidFluidSolver` (Cases 2, 3, free surface) with static `matrix/rhs/solve/dispersion/unknowns`; `Interface` dispatches | requested names; M built once per (k, ω) and solved for all φ with `pagemldivide` |
| `findPoles` → real poles | `[poles, leaky] = findPoles(...)` also returns **near-real leaky poles** (real-axis Newton) | slow fluid (α_f < c_R): leaky Rayleigh peak of width ∝ ρ_f α_f is otherwise missed by the quadrature *and* by the balance check |
| `windowNodes(breaks, n)` | `windowNodes(breaks, n, peaks)` grades nodes at a ± b·10^j around leaky peaks | resolves Lorentzians of any width (tested down to b ≈ 1e-13) |
| `EP.pct` | adds `.leaky`, `.trapped` (= real IW + quasi-trapped leaky) | Rayleigh limit of a slow fluid |
| `synthetics` only | + `wavefield/freq2time.m` (complex-ω → time, `fft` kernel, `exp(eps t)`) | testable separately (Ricker round trip) |
| DWN Bessel sum | + Euler–Maclaurin end correction at k = 0 (n = 0 harmonic) | trapezoid on [0, ∞) of k·u₀·J₀ is O(Δk²) otherwise (1e-4 errors observed) |
| `sourcePower(iface, src, poles, ω, opts)` | `sourcePower(iface, src, ω, W0, Wi, leaky, opts)` | reuses W0 and pole terms already computed |
| `pattern` option | not implemented | not needed yet |

Run (from the repository root): `runtests('tests/lambfk/testLamb.m'); run('tests/lambfk/test_rayleigh_energy_limit.m')`
