function K = bornKernelsRayleigh(M, z, opts)
%BORNKERNELSRAYLEIGH  Born coupling kernels for delta-rho, delta-lambda, delta-mu.
%
%   K = bornKernelsRayleigh(M, z) returns, at scatterer depth z, the
%   angle-averaged Born coupling kernels between the Rayleigh mode of M and a
%   bulk P or S wave, computed TWO INDEPENDENT WAYS:
%
%   K.xu.*     path A -- the kernels Xu et al. (2021) publish in Appendix B,
%              transcribed directly (B.1 for delta-lambda, B.3-B.5 and
%              B.12-B.14 for delta-rho and delta-mu).
%
%   K.born.*   path B -- computed from the general Born matrix element
%
%                M_ba = INT [ w^2 drho u_b* . u_a
%                             - dlambda (div u_b*)(div u_a)
%                             - 2 dmu   eps_b* : eps_a ] dV
%
%              with u_a the Rayleigh eigenfunction and u_b a bulk plane wave,
%              angle-averaged over the bulk direction numerically. Nothing
%              from Xu's appendix enters this path; only the eigenfunctions.
%
%   The two paths share no algebra, so agreement is a real cross-check of
%   both the transcription and the derivation.
%
%   Sign convention: Xu's z points DOWN. The Born path is built in the same
%   frame, so the cross terms r1*dr2 and dr1*r2 come out with Xu's signs. If
%   the convention is flipped those two terms change sign and the squares do
%   not -- which is exactly what this comparison is able to detect.
%
%   opts.nDir  directions in the angular average (default 40000)
%
%   TYPO IN Xu et al. (2021) B.11, B.14 and B.22. The two paths agree on
%   every kernel except mu_Sh, and the evidence says the published B.11 is
%   the one at fault:
%
%     printed B.11 : r1^2 sin^2(psi) + r2^2 cos^2(psi)
%                    + (dr1/kR)^2 - 2 (dr1/kR) r2 cos(psi)
%     Born route   : r1^2 sin^2(psi) + cos^2(psi) * [ dr1/kR - r2 ]^2
%
%   Three independent reasons, all checked against the typeset paper (not
%   the OCR, which transcribed B.11 correctly):
%
%   1. The same Born machinery reproduces B.9 and B.10 EXACTLY, angle by
%      angle -- twelve terms carrying sin^4, cos^4, sin^2cos^2, sin^2(2psi)
%      and cos^2(2psi). Ratios are 1.00000000 at psi = 10, 35, 60, 85, 120
%      degrees. Only B.11 fails.
%   2. B.9 and B.10 are even under psi -> pi-psi (0.7422911 and 0.3547244 at
%      both 60 and 120 degrees). Printed B.11 is not (0.1179 vs 0.8779),
%      because cos(psi) appears to the first power. The three kernels
%      describe the same geometry and must share that symmetry.
%   3. B.9 and B.10 both end in [angular factor] * (dr1/kR - r2)^2. Printed
%      B.11 breaks that pattern; the Born form restores it.
%
%   So two angular factors are missing from the printed B.11: cos^2(psi) on
%   the (dr1/kR)^2 term, and cos(psi) should be cos^2(psi) on the cross
%   term. This propagates: B.14 and B.22 are the angular average of B.11, so
%   their brackets change from
%       3 (dr1/kR)^2            to    (dr1/kR)^2 - 2 (dr1/kR) r2
%   K.xu.mu_Sh below is the PRINTED form, kept so the discrepancy stays
%   visible; K.xu.mu_Sh_fixed is the corrected one, which equals the Born
%   route. Nothing in this repository consumed mu_Sh, so no result changes.
%
%   See also RAYLEIGHMODEL, TESTBORNKERNELS.
%
%   References:
%     - Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Monte Carlo simulations of coupled
%       body- and Rayleigh-wave multiple scattering in elastic media, Geophys. J. Int. 228,
%       1213-1236. doi:10.1093/gji/ggab370

if nargin < 3, opts = struct; end
if ~isfield(opts,'nDir') || isempty(opts.nDir), opts.nDir = 40000; end

kR = M.kR;  cR = M.cR;
z = z(:);

% --- eigenfunctions and their depth derivatives (4th-order, then checked
%     against the analytic form in testBornKernels) ---------------------
h  = 1e-4/kR;
d1 = @(f,zz) (-f(zz+2*h) + 8*f(zz+h) - 8*f(zz-h) + f(zz-2*h))/(12*h);
r1  = M.r1(z);      r2  = M.r2(z);
dr1 = d1(M.r1,z);   dr2 = d1(M.r2,z);

% =====================================================================
% PATH A: Xu et al., Appendix B
% =====================================================================
K.xu.rho_P  = r1.^2 + r2.^2;                       % B.6 kernel (angle-averaged)
K.xu.rho_Sv = r1.^2 + 4*r2.^2;                     % B.7
K.xu.rho_Sh = r1.^2;                               % B.8
K.xu.lam_P  = (r1 + dr2/kR).^2;                    % B.1
K.xu.lam_Sv = zeros(size(z));                      % B.2
K.xu.lam_Sh = zeros(size(z));                      % B.2
K.xu.mu_P   = 3*r1.^2 + 3*(dr2/kR).^2 + 2*r1.*dr2/kR ...
              + (dr1/kR).^2 + r2.^2 - 2*(dr1/kR).*r2;   % B.12 / B.20 bracket
K.xu.mu_Sv  = r1.^2/20 + (2/15)*(dr2/kR).^2 - (2/15)*r1.*dr2/kR ...
              + (7/120)*r2.^2 + (7/120)*(dr1/kR).^2 ...
              - (7/60)*(dr1/kR).*r2;                    % B.13 / B.21 bracket
K.xu.mu_Sh  = r1.^2/12 + r2.^2/24 + (1/8)*(dr1/kR).^2;  % B.14 / B.22 AS PRINTED
% corrected: angular average of r1^2 sin^2 + cos^2*(dr1/kR - r2)^2, scaled
% to the same 1/12 leading coefficient as the printed version
K.xu.mu_Sh_fixed = r1.^2/12 + (1/24)*((dr1/kR) - r2).^2;

% =====================================================================
% PATH B: general Born matrix element, angle-averaged numerically
% =====================================================================
% Rayleigh displacement in Xu's frame (z down), travelling along +x. The
% relative sign of the vertical component is NOT free: it is fixed by the
% traction-free surface, and checked here rather than assumed. With
% RayleighModel's r1, r2,
%
%   shear   r1'(0) - kR r2(0) = 0     (relative residual 1.9e-10)
%   normal  lam(-kR r1 - r2') - 2 mu r2' = 0   (8.8e-11)
%
% both hold for
%
%   u = ( i r1, 0, -r2 ) exp(i kR x)
%
% and neither holds for +r2. Getting this wrong is invisible in delta-rho --
% only r1^2 and r2^2 appear there -- but it flips the cross terms r1*dr2 and
% dr1*r2 that delta-lambda and delta-mu depend on. testBornKernels asserts
% the two boundary conditions so the convention cannot drift.
exx = -kR*r1;   ezz = -dr2;   exz = 0.5i*(dr1 - kR*r2);
divR = -kR*r1 - dr2;                                % = -kR (r1 + dr2/kR)

% quasi-uniform directions on the sphere (Fibonacci), for the angular average
nD = opts.nDir;
ii = (0:nD-1).' + 0.5;
nz_ = 1 - 2*ii/nD;
rr = sqrt(max(0,1-nz_.^2));
ga = pi*(1+sqrt(5));
nx_ = rr.*cos(ga*ii);   ny_ = rr.*sin(ga*ii);

nP = numel(z);
K.born.rho_P  = zeros(nP,1);  K.born.rho_Sv = zeros(nP,1);
K.born.rho_Sh = zeros(nP,1);  K.born.lam_P  = zeros(nP,1);
K.born.lam_Sv = zeros(nP,1);  K.born.lam_Sh = zeros(nP,1);
K.born.mu_P   = zeros(nP,1);  K.born.mu_Sv  = zeros(nP,1);
K.born.mu_Sh  = zeros(nP,1);

for q = 1:nP
    uR = [1i*r1(q); 0; -r2(q)];                    % Rayleigh displacement
    E  = [exx(q) 0 exz(q); 0 0 0; exz(q) 0 ezz(q)];% Rayleigh strain

    % ---- P wave: e_b = n, div = i kP, eps = i kP n(x)n ----------------
    dotP = nx_*uR(1) + ny_*uR(2) + nz_*uR(3);      % e_P . u_R
    % The 3, 6 and 2 below are Xu's per-channel normalisations, not fudge:
    % the three direction-averaged overlaps sum to |u|^2 by completeness of
    % the polarization triad, (r1^2+r2^2)/3 + (r1^2+4r2^2)/6 + r1^2/2
    % = r1^2 + r2^2, and B.6-B.8 quote them each scaled to unit leading term.
    K.born.rho_P(q) = 3*mean(abs(dotP).^2);

    % eps_P* : eps_R  (drop the common i*kP, restored by the calibration)
    cc = nx_.^2*E(1,1) + nz_.^2*E(3,3) + 2*nx_.*nz_*E(1,3);
    K.born.mu_P(q) = 15/kR^2*mean(abs(cc).^2);     % 15/kR^2 : Xu's B.12 scale

    K.born.lam_P(q) = mean(abs(divR(q)).^2)/kR^2;  % (div u_R)^2 / kR^2

    % ---- S waves: two transverse polarizations per direction ----------
    % Build an orthonormal transverse pair (SV in the n-z plane, SH normal).
    shx = -ny_;  shy = nx_;  shn = hypot(shx,shy);
    deg = shn < 1e-12;  shx(deg) = 1; shy(deg) = 0; shn(deg) = 1;
    shx = shx./shn;  shy = shy./shn;
    svx = shy.*nz_;            svy = -shx.*nz_;
    svz = shx.*ny_ - shy.*nx_;

    dSv = svx*uR(1) + svy*uR(2) + svz*uR(3);
    dSh = shx*uR(1) + shy*uR(2);
    K.born.rho_Sv(q) = 6*mean(abs(dSv).^2);
    K.born.rho_Sh(q) = 2*mean(abs(dSh).^2);

    % div of a shear wave vanishes identically -> no delta-lambda coupling
    K.born.lam_Sv(q) = mean(abs(nx_.*svx + ny_.*svy + nz_.*svz).^2)*abs(divR(q))^2/kR^2;
    K.born.lam_Sh(q) = mean(abs(nx_.*shx + ny_.*shy).^2)*abs(divR(q))^2/kR^2;

    % delta-mu into a shear wave: eps_b = (i k/2)(n(x)e + e(x)n), so
    % eps_b* : eps_R = -i k * n' * eps_R * e. The constants below are Xu's
    % per-channel normalisations, the same role the 3/6/2 play for delta-rho.
    nEn = @(ex,ey,ez) (nx_.*E(1,1) + nz_.*E(3,1)).*ex ...
                    + (nx_.*E(1,2) + nz_.*E(3,2)).*ey ...
                    + (nx_.*E(1,3) + nz_.*E(3,3)).*ez;
    K.born.mu_Sv(q) = 1/kR^2*mean(abs(nEn(svx,svy,svz)).^2);
    K.born.mu_Sh(q) = 1/kR^2*mean(abs(nEn(shx,shy,0*shx)).^2);
end
end
