# MathDev — Extended Lamb's Problem: f–k Formulation & Energy Partition

## Goal
MATLAB code for the elastodynamic Green's function of a two-half-space interface in the
frequency–wavenumber (f–k) domain, computing the **energy partition (%)** between
body waves (P/S reflected/transmitted, fluid P) and interface waves (Rayleigh/Scholte/Stoneley).

**Phase 1 (now):** full mathematical formulation in LaTeX (`LambFK_Formulation.tex`). No MATLAB yet.
**Phase 2 (later):** MATLAB implementation, after the formulation is approved.

## Required content (from the task)
- [x] Boundary-condition systems `M x = S` in the f–k domain for:
  - [x] Case 1: solid–solid, welded contact, source in solid 1 (4×4)
  - [x] Case 2: solid–fluid, source in the solid (3×3)
  - [x] Case 3: solid–fluid, source in the fluid (3×3)
- [x] Solid–fluid 3×3 rows: continuity of u_z, continuity of σ_zz (= −p), σ_xz = 0
- [x] Energy: % of radiated energy in the interface wave (pole residue) vs body waves
      (real-axis / branch-cut integral)
- [x] Degeneration of solid–fluid to free-surface (Rayleigh) when ρ_f, α_f → 0
- [x] LaTeX must compile cleanly

## Conventions (fixed — keep consistent in MATLAB)
- 3D target: point force (any direction) or arbitrary moment tensor. Transform (kx,ky) → polar (k,φ);
  per (k,φ) the problem = 2D P-SV in plane (l̂,ẑ) with the SAME M1/M2 + SH (2×2 solid–solid, C1=C0 fluid/free).
  2D line-source sections kept as building blocks.
- Time dependence e^{−iωt}; transform f̃(k) = ∫ f(x) e^{−ikx} dx, inverse (1/2π)∫ … e^{ikx} dk.
- Interface at z = 0, z positive down. Solid (medium 1, source side in Cases 1–2) occupies z > 0;
  medium 2 / fluid occupies z < 0.
- Potentials: u = ∇φ + ∇×(ψ ŷ) → u_x = ∂xφ − ∂zψ, u_z = ∂zφ + ∂xψ.
- Vertical wavenumbers ν_m = sqrt(k² − k_m²), k_m = ω/c_m, Re ν ≥ 0, Im ν ≤ 0.
  MATLAB-safe: `nu = -1i*sqrt(km^2 - k.^2)`.
- γ = 2k² − k_β²; Rayleigh function R(k) = γ² − 4k²ν_α ν_β.
- Same M for a given interface regardless of source side — only S changes (Cases 2 & 3 share M).

## Added requirements (2026-09-23, round 2)
- [x] 3D domain; point force any direction; arbitrary moment tensor (tex §7)
- [x] MATLAB architecture: classdef/structs, energy module (R/T coeffs + residues), signatures+docstrings only
      → `MATLAB_Architecture.md`
- DWN + inverse FFT with ω_c = ω + iε for wavefields; energy partition evaluated at REAL ω (decided).

## Energy method (decided)
- Time-averaged flux ⟨P_j⟩ = (ω/2) Im(σ_ij u_i*).
- Body waves: vertical flux through far horizontal planes via Parseval;
  per-unit-amplitude flux E_m(k) = ½ ρ ω³ Re(i ν_m) (nonzero only for |k| < k_m).
  Interface poles lie outside the propagating windows → integrals are regular.
- Interface wave: residue at real pole k_p (Res M⁻¹ = r lᵀ/(lᵀ M′ r)), horizontal flux
  integrated over depth; both sides x → ±∞.
- Leaky poles (leaky Rayleigh, leaky Stoneley) are NOT trapped energy → counted as body waves.
- Percentages: η = W_part / W_total.

## Validation checks (to carry into MATLAB)
1. Energy balance: W_source (work done by source, Sokhotski–Plemelj at the pole)
   = W_body + W_interface.
2. Plane-wave check: energy R/T coefficients sum to 1 for propagating incidence.
3. ρ_f → 0 recovers free-surface Rayleigh pole and 2×2 free-surface system.
4. h = 0 (surface source): partition independent of frequency.
5. 3D extension benchmark: Miller & Pursey (1955), Poisson solid, vertical surface force:
   ≈ 67% Rayleigh, 26% S, 7% P.

## Status / log
- 2026-09-23: folder was empty; formulation started; this file created.
- 2026-09-23: `LambFK_Formulation.tex` written, compiles (MacTeX; `export PATH=/Library/TeX/texbin:$PATH`,
  `latexmk -pdf -outdir=out`). Energy formulas verified numerically (scratch Python, free-surface limit,
  vertical force): W_src = W_body + W_IW to ~1e-9 relative error.
  Reference partition (Poisson solid, 2D line vertical force, h = 0): P 21.94 %, S 29.13 %, Rayleigh 48.93 %
  — use as MATLAB regression target.
- Gotcha found: residue by finite difference needs (k−k_p) with its sign at k_p−δ; a sign slip
  gives zero residue AND still passes the energy balance (both sides lose the pole term) → check η_IW > 0.
- Next: approve formulation → Phase 2 MATLAB (Case 2 free-surface limit first, then Cases 2/3, then Case 1).
- 2026-09-23 (round 2): 3D section added to tex (Weyl kernel same as 2D; direct field from Green's tensor
  via κ vectors; Π/Σ/H amplitudes; SH systems; W_body (eq:Wbody3), W_IW = k_p/2π ∫ W2D+ dφ, source work
  with M*:∇u). Verified numerically: whole-space force power vs analytic; S/P ratios 2(α/β)³ and 1.5(α/β)⁵;
  Miller–Pursey 6.88/25.75/67.36 %; energy balance ~1e-8 for general M and oblique F at k_βh = 0.5.
- Gotcha: body integrals have 1/sqrt singularities at INTERIOR branch points (k_α inside S window,
  grazing-P→SV/head waves). Must split windows at all branch points + sin² substitution; midpoint rule
  gave erratic 1e-3 errors. Earlier 2D claim "no singularities" corrected in tex.
- Gotcha: with e^{-iωt}, synthesis uses MATLAB `fft` (not `ifft`) and ω_c = ω + iε (upward shift).
- `MATLAB_Architecture.md` written. Next: approve → implement in order of its §7.
- 2026-09-24: MATLAB package implemented (`matlab/+lamb`, see MATLAB_Architecture.md §8). Tests 7/7 pass:
  whole-space power & S/P ratios (1e-8), Miller–Pursey, energy balance for all 3 cases + free surface
  (force, MT, explosion, fluid source; ~1e-10), R/T unitarity (1e-10), weak-fluid limit, DWN vs
  analytic Green's tensor at complex ω (force 1e-6, MT 1e-5), Ricker complex-ω round trip (1e-8).
  End-to-end Lamb synthetics: Rayleigh arrival at r/c_R, √r spreading.
- PHYSICS CORRECTION: tex had the leaky condition reversed. Leaky Rayleigh occurs for SLOW fluid
  (α_f < c_R, e.g. water/rock, α_f → 0); fast fluid gives a real Scholte pole → Rayleigh. For α_f → 0
  the partition limit is non-uniform (r → ∞: leaky energy ends in fluid P). Fixed in tex §4.
- Gotcha: a near-real leaky pole (b ~ 1e-13) is missed by BOTH body quadrature and source-work
  integral → energy balance closes at 1e-16 while the partition is wrong (IW 0 %). Balance check
  alone cannot detect missing poles; always also check against a known limit.
- Gotcha: cylindrical DWN trapezoid needs the k=0 Euler–Maclaurin correction (+Δk²/12·u₀(0)).
- `test_rayleigh_energy_limit.m`: ρ_f = α_f = 1e-6 → trapped 67.36 %, P 6.88 %, S 25.75 % ✓.
- 2026-09-24: Zoeppritz validation DONE.
  (a) plane wave: lamb.energy.rtCoefficients == MaterialClass.ZoeppritzSolid for P, SV, SH incidence,
      0.5–89.5° (incl. post-critical), max diff ~1e-14 (shale/sandstone 2500/1100/2200 | 3000/1700/2300).
  (b) full waveform `tests/test_zoeppritz_point_source.m`: explosion at H = 100 λ, θ = 30°, Ricker 25 Hz,
      receivers on specular rays, windowed ρc∫|v|²dt, exact plane-interface ray-tube spreading
      X|dX/di|cos(i_out)/(R² sin i). Rpp 0.1246 / Rps 1.6041 / Tpp 94.765 / Tps 3.506 % vs Zoeppritz:
      all within 0.01 % relative; tolerance set to 0.5 %. Runtime ~80 s.
  Added `dwn`/`synthetics` option 'Direct' (scattered field only) and wavelet-band frequency skipping.
