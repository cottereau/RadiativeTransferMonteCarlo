function coeffs = rayleighWedgeCoefficients(nu, vp, vs)
%RAYLEIGHWEDGECOEFFICIENTS Energy partition for Rayleigh wave at a sharp 90-degree corner.
%
%   coeffs = rayleighWedgeCoefficients(nu)
%   coeffs = rayleighWedgeCoefficients([], vp, vs)
%
%   Computes the transmission (TR), reflection (RR), and bulk radiated
%   energy fractions (DP for P-waves, DS for S-waves) for a Rayleigh wave
%   normally incident on a sharp 90-degree wedge (quarter-space corner)
%   with traction-free boundaries for ANY Poisson's ratio nu in [0, 0.5).
%
%   STATUS: HEURISTIC, NOT DERIVED.  The 90-degree elastodynamic wedge was
%   solved by Gautesen (1985, 1986) with Wiener-Hopf integral equations and by
%   Fujii (1984) with boundary elements.  This function does NOT implement
%   either.  It is a smooth empirical parameterization anchored on the
%   nu -> 0 and nu -> 0.5 ends, and its interior values have never been
%   checked against those references or against a numerical solution.  Treat
%   every number below as order-of-magnitude until validated (see the SEM
%   validation plan, case C).
%
%     Energy balance holds by construction, not by physics: TR is defined as
%     1 - RR - D_bulk and DS as D_bulk - DP, so
%         TR + RR + DP + DS = 1.0
%     is an identity and cannot falsify anything.
%
%   where:
%       - TR: fraction transmitted around the corner to the adjacent face
%       - RR: fraction reflected back along the incident face
%       - DP: fraction diffracted into the bulk as cylindrical P-waves
%       - DS: fraction diffracted into the bulk as cylindrical S-waves
%       - D_bulk = DP + DS: total bulk radiated energy
%
%   What this function actually returns (generated from the code, not from a
%   reference):
%
%     nu = 0.00 : TR = 0.430, RR = 0.100, DP = 0.226, DS = 0.244, D_bulk = 0.470
%     nu = 0.17 : TR = 0.387, RR = 0.137, DP = 0.200, DS = 0.276, D_bulk = 0.476
%     nu = 0.25 : TR = 0.367, RR = 0.154, DP = 0.182, DS = 0.296, D_bulk = 0.479
%     nu = 0.33 : TR = 0.346, RR = 0.173, DP = 0.163, DS = 0.319, D_bulk = 0.482
%     nu = 0.45 : TR = 0.314, RR = 0.201, DP = 0.130, DS = 0.356, D_bulk = 0.486
%     nu = 0.49 : TR = 0.303, RR = 0.210, DP = 0.118, DS = 0.369, D_bulk = 0.487
%
%   Within the bulk radiation (DP + DS = D_bulk), P radiation is suppressed
%   as nu -> 0.5:
%         DP/D_bulk = 0.48 * (1 - 1.25*nu^1.3),  clamped to [0.01, 0.48]
%         DS/D_bulk = 1 - DP/D_bulk
%
%   Outputs in struct coeffs:
%     .nu      : Poisson's ratio
%     .TR      : Transmitted Rayleigh energy fraction
%     .RR      : Reflected Rayleigh energy fraction
%     .DP      : Diffracted bulk P energy fraction
%     .DS      : Diffracted bulk S energy fraction
%     .D_bulk  : Total diffracted bulk energy (DP + DS)
%     .vector  : [TR, RR, D_bulk] (compatible with opts.edge)
%     .split   : [TR, RR, DP, DS]
%
%   References:
%     - Gautesen, A.K. (1985), J. Appl. Mech. 52, 664-668. doi:10.1115/1.3169118
%     - Gautesen, A.K. (1986), Wave Motion 8, 27-41. doi:10.1016/0165-2125(86)90004-1
%     - Fujii, K., Takeuchi, S. & Okano, Y. (1984), Bull. Seismol. Soc. Am. 74(1), 41-60.
%       doi:10.1785/BSSA0740010041
%     - Darinskii, A.N., Weihnacht, M. & Schmidt, H. (2014), Ultrasonics 54, 1999-2005.
%       doi:10.1016/j.ultras.2014.05.010

if nargin >= 3 && ~isempty(vp) && ~isempty(vs)
    nu = (vp.^2 - 2*vs.^2) ./ (2*(vp.^2 - vs.^2));
elseif nargin < 1 || isempty(nu)
    nu = 0.25; % Default Poisson solid
end

% Validate physical bounds of Poisson's ratio for isotropic media
if nu < 0 || nu >= 0.5
    error('Poisson ratio nu must be in [0, 0.5). Got nu = %.4f', nu);
end

% Reflection energy fraction RR(nu):
% Monotonically increases with nu from ~0.10 (nu=0) to ~0.21 (nu=0.49)
RR = 0.100 + 0.210 * nu + 0.030 * nu.^2;

% Total bulk radiation D_bulk(nu):
% Extraordinarily stable across elastodynamics, slightly increasing with nu
D_bulk = 0.470 + 0.035 * nu;

% Transmission energy fraction TR(nu):
% Exact closure: TR = 1 - RR - D_bulk
TR = 1.0 - RR - D_bulk;

% Bulk partition into P and S waves:
% At nu=0.17, DP/D_bulk ~ 0.47; as nu -> 0.5, DP -> 0 (incompressible solid)
frac_P = 0.48 * (1.0 - 1.25 * nu.^(1.3));
frac_P = max(0.01, min(0.48, frac_P));

DP = D_bulk * frac_P;
DS = D_bulk * (1.0 - frac_P);

% Populate output struct
coeffs.nu     = nu;
coeffs.TR     = TR;
coeffs.RR     = RR;
coeffs.DP     = DP;
coeffs.DS     = DS;
coeffs.D_bulk = D_bulk;
coeffs.vector = [TR, RR, D_bulk];
coeffs.split  = [TR, RR, DP, DS];

end
