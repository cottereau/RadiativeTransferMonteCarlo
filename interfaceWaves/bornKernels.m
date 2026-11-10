function K = bornKernels(ux, uz, dux, duz, k, nDir)
%BORNKERNELS  Born coupling kernels of an interface mode, frame-independent.
%
%   K = bornKernels(ux, uz, dux, duz, k) returns the angle-averaged Born
%   kernels coupling an interface mode to bulk P, Sv and SH waves, for
%   perturbations in density, in lambda and in mu.
%
%   Inputs are the mode's own complex displacement components and their
%   depth derivatives at the scatterer, in WHATEVER frame the model uses,
%   with the mode travelling along +x:
%
%       u = ( ux(z), 0, uz(z) ) exp(i k x)
%
%   Nothing here assumes Xu's (r1, r2) convention or which way z points.
%   That matters: the relative sign of the horizontal and vertical
%   components is invisible in delta-rho but flips the cross terms that
%   delta-lambda and delta-mu are built from, and getting it wrong once
%   produced factor-1600 errors (see testBornKernels). Each model's own sign
%   is already fixed by the interface conditions it satisfies, so feeding
%   the model's own handles in is safer than translating conventions.
%
%   Normalisation follows Xu et al. (2021) Appendix B channel by channel, so
%   for the Rayleigh eigenfunctions these reproduce B.1-B.14 exactly -- with
%   the single exception of B.11, whose printed form carries two wrong
%   angular factors (bornKernelsRayleigh.m documents the evidence). K.mu_Sh
%   here is the corrected kernel.
%
%   K.rho_P, K.rho_Sv, K.rho_Sh, K.lam_P, K.mu_P, K.mu_Sv, K.mu_Sh
%   (lam_Sv and lam_Sh vanish identically: a shear wave has no divergence.)
%
%   See also BORNKERNELSRAYLEIGH, SCHOLTEMODEL, STONELEYMODEL.
%
%   References:
%     - Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Monte Carlo simulations of coupled
%       body- and Rayleigh-wave multiple scattering in elastic media, Geophys. J. Int. 228,
%       1213-1236. doi:10.1093/gji/ggab370

if nargin < 6 || isempty(nDir), nDir = 40000; end
ux = ux(:); uz = uz(:); dux = dux(:); duz = duz(:);
n = numel(ux);

% strain of the interface mode
exx = 1i*k*ux;
ezz = duz;
exz = 0.5*(1i*k*uz + dux);
divu = exx + ezz;

% quasi-uniform directions on the sphere (Fibonacci)
ii = (0:nDir-1).' + 0.5;
nz_ = 1 - 2*ii/nDir;
rr  = sqrt(max(0,1-nz_.^2));
ga  = pi*(1+sqrt(5));
nx_ = rr.*cos(ga*ii);  ny_ = rr.*sin(ga*ii);

% transverse pair per direction: SH horizontal, SV completing the triad
shx = -ny_;  shy = nx_;  shn = hypot(shx,shy);
deg = shn < 1e-12;  shx(deg) = 1; shy(deg) = 0; shn(deg) = 1;
shx = shx./shn;  shy = shy./shn;
svx = shy.*nz_;  svy = -shx.*nz_;  svz = shx.*ny_ - shy.*nx_;

K.rho_P = zeros(n,1); K.rho_Sv = zeros(n,1); K.rho_Sh = zeros(n,1);
K.lam_P = zeros(n,1); K.mu_P  = zeros(n,1);
K.mu_Sv = zeros(n,1); K.mu_Sh = zeros(n,1);

for q = 1:n
    U = [ux(q); 0; uz(q)];
    E = [exx(q) 0 exz(q); 0 0 0; exz(q) 0 ezz(q)];

    % delta-rho: overlap of the bulk polarization with the mode displacement
    K.rho_P(q)  = 3*mean(abs(nx_*U(1) + nz_*U(3)).^2);
    K.rho_Sv(q) = 6*mean(abs(svx*U(1) + svz*U(3)).^2);
    K.rho_Sh(q) = 2*mean(abs(shx*U(1)).^2);

    % delta-lambda: only the divergence, so only P couples
    K.lam_P(q) = abs(divu(q))^2/k^2;

    % delta-mu: eps_b* : eps_mode = -i k_b * n' * eps_mode * e
    nEn = @(ex,ey,ez) (nx_*E(1,1) + nz_*E(3,1)).*ex + (nx_*E(1,3) + nz_*E(3,3)).*ez;
    K.mu_P(q)  = 15/k^2*mean(abs(nEn(nx_, ny_, nz_)).^2);
    K.mu_Sv(q) =  1/k^2*mean(abs(nEn(svx, svy, svz)).^2);
    K.mu_Sh(q) =  1/k^2*mean(abs(nEn(shx, shy, 0*shx)).^2);
end
end
