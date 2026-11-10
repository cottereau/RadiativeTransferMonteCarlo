function M = RayleighModel(alpha0,beta0,freq,V,drho,n,zmax)
%RAYLEIGHMODEL  Coupled body/Rayleigh scattering model (Xu et al., 2021).
%
%   M = RayleighModel(alpha0,beta0,freq,V,drho,n,zmax) builds everything
%   propagateParticleRayleigh needs for a homogeneous elastic halfspace
%   seeded with point density scatterers (delta-rho Born regime):
%   Rayleigh velocity and eigenfunctions (Appendix A), angle-averaged
%   cross sections (Appendix B) and the injection-depth tables (Appendix D)
%   of
%
%     Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Monte Carlo simulations
%     of coupled body- and Rayleigh-wave multiple scattering in elastic
%     media, Geophys. J. Int. 228, 1213-1236. doi:10.1093/gji/ggab370
%
%   Modes are numbered   1 = P,  2 = Sv,  3 = Sh,  4 = Rayleigh.
%   Depth z is positive downwards, the free surface is z = 0.
%
%   Units are arbitrary but must be consistent. With alpha0, beta0 in km/s,
%   V in km^3 and n in km^-3 the Xu et al. benchmark returns l_R = 35.8 km.
%   (The "0.5 m^3 / 1 m^-3" printed in their Section 4 is a typography
%   slip: the scattering rate goes as n*V^2 and is not scale invariant.)
%
%   See also PROPAGATEPARTICLERAYLEIGH, VALIDATIONXU2021, TESTRAYLEIGHMODEL.

if nargin < 7 || isempty(zmax), zmax = 10; end

w = 2*pi*freq;

% ---- Rayleigh phase velocity: the real root of the secular equation ------
secular = @(c) (2-(c/beta0).^2).^2 ...
    - 4*sqrt(max(0,1-(c/alpha0).^2)).*sqrt(max(0,1-(c/beta0).^2));
cR = fzero(secular,[0.5*beta0 0.9999*beta0]);
UR = cR;                    % no dispersion in a homogeneous halfspace

kR = w/cR;  kP = w/alpha0;  kS = w/beta0;
ea = sqrt(1-(cR/alpha0)^2);
eb = sqrt(1-(cR/beta0)^2);

% ---- displacement eigenfunctions (Xu et al. A.3-A.4) --------------------
r1 = @(z) -kR*exp(-kR*ea*z) + kR*((2*beta0^2-cR^2)/(2*beta0^2))*exp(-kR*eb*z);
r2 = @(z) -kR*ea*exp(-kR*ea*z) + kR*(2*ea/(2-(cR/beta0)^2))*exp(-kR*eb*z);

zinf = 50/(kR*min(ea,eb));  % 50 e-folding depths stands in for infinity
I    = @(f) integral(f,0,zinf,'RelTol',1e-11,'AbsTol',0);

Inorm  = I(@(z) r1(z).^2 + r2(z).^2);          % energy normalisation
Ir1    = I(@(z) r1(z).^2);
Ir14r2 = I(@(z) r1(z).^2 + 4*r2(z).^2);
Ir4    = I(@(z) r1(z).^4 + 2*r2(z).^4);
Ir2    = Inorm - Ir1;

e2 = V^2*drho^2;            % V^2 (delta_rho/rho0)^2, common to every cross section

% ---- body-to-body cross sections, averaged over angle (B.28-B.39) -------
% ponytail: totals only, so the walk scatters isotropically. Equipartition
% follows from the totals alone (detailed balance), so the benchmark is
% unaffected; add the psi0-dependent phase functions B.29-B.33 if transport
% mean free paths or angular distributions are ever needed.
sPP   = kP^4/(12*pi)*e2;
sPS   = beta0*kS^4/(12*pi*alpha0)*e2;         % identical for Sv and Sh
sSP   = alpha0*kP^4/(12*pi*beta0)*e2;         % identical from Sv and Sh
sSSl  = kS^4/(8*pi)*e2;                       % Sv->Sv and Sh->Sh
sSSx  = kS^4/(24*pi)*e2;                      % Sv->Sh and Sh->Sv
M.sigBB = [sPP  sPS   sPS ; ...
           sSP  sSSl  sSSx; ...
           sSP  sSSx  sSSl];

% ---- Rayleigh-to-body (B.17-B.19) and Rayleigh-to-Rayleigh (B.24) ------
M.sigRB = [ w*kP^3/(12*pi*UR)*e2 , ...
            w*kS^3/(24*pi*UR)*(Ir14r2/Inorm)*e2 , ...
            w*kS^3/( 8*pi*UR)*(Ir1   /Inorm)*e2 ];
M.sigRR = w^3/(8*cR*UR^2*Inorm^2)*Ir4*e2;

% ---- body-to-Rayleigh, depth and angle dependent (B.3-B.5) -------------
% sigma_PR  = cP*[ r1^2 sin^2(psi0) + 2 r2^2 cos^2(psi0) ]
% sigma_SvR = cS*[ r1^2 cos^2(psi0) + 2 r2^2 sin^2(psi0) ]
% sigma_ShR = cS*  r1^2
% with psi0 the incidence angle from the vertical, so that cos(psi0) is the
% vertical cosine of the particle and the two depth profiles r1^2 and r2^2
% are the only tables needed. Keeping the angle instead of the B.6-B.8
% average is what makes the coupling reciprocal for an anisotropic
% body-wave field, which is exactly the situation near the free surface.
dz = zmax/2000;
zg = (0:dz:zmax).';
nz = numel(zg);
M.cP = w^3/(8*alpha0*cR*UR*Inorm)*e2;
M.cS = w^3/(8*beta0 *cR*UR*Inorm)*e2;

% vertical optical depths n*int_0^z r1^2 dz' and n*int_0^z r2^2 dz', so a
% straight step integrates the depth dependence exactly
M.G1 = cumtrapz(zg,n*r1(zg).^2);
M.G2 = cumtrapz(zg,n*r2(zg).^2);

% angle-averaged version (B.6-B.8), kept for reciprocity checks and plots
sigBR = [ 2*M.cP/3*(r1(zg).^2 + r2(zg).^2) , ...
            M.cS/3*(r1(zg).^2 + 4*r2(zg).^2) , ...
            M.cS  * r1(zg).^2 ];
M.FBR = cumtrapz(zg,n*sigBR);

% ---- Rayleigh-to-body emission kernel (Appendices D and E) -------------
% By reciprocity the outgoing (depth, colatitude) pair follows the same
% kernel, which splits into two separable components: one with z ~ r1^2 and
% one with z ~ r2^2, each with its own law for mu = cos(psi1).
%   mu-law 0 : uniform      1 : proportional to 1-mu^2      2 : to mu^2
M.emitLaw = [1 2; 2 1; 0 0];              % rows: P, Sv, Sh
wA = [Ir1*4/3 2*Ir2*2/3;                  % P
      Ir1*2/3 2*Ir2*4/3;                  % Sv
      Ir1       0      ];                 % Sh
M.emitP1 = wA(:,1)./sum(wA,2);            % probability of the r1^2 component

zf  = linspace(0,zmax,20001).';
M.zInv = cell(1,2);
pdfz = [r1(zf).^2 , r2(zf).^2];
for c = 1:2
    cdf  = cumtrapz(zf,pdfz(:,c));
    cdf  = cdf/cdf(end);
    keep = [true; diff(cdf) > 0];   % r1 vanishes near 0.5 km: flat CDF there
    M.zInv{c} = griddedInterpolant(cdf(keep),zf(keep),'linear','nearest');
end

% inverse CDF of mu ~ 1-mu^2 on [-1,1] (the mu^2 law is inverted in closed
% form in the propagator, and the uniform one needs nothing)
mg   = linspace(-1,1,20001).';
cdfm = (3*mg - mg.^3 + 2)/4;
keep = [true; diff(cdfm) > 0];
M.muInv = griddedInterpolant(cdfm(keep),mg(keep),'linear','nearest');

% ---- azimuth of the conversion, B.3-B.5 and B.9-B.11 -------------------
% The same two-component split that carries the depth and the polar angle
% also carries the azimuth. Reading B.3-B.5 with the incidence plane as x-z
% and the interface wave leaving at azimuth phi,
%
%   P  : r1^2 sin^2(psi) cos^2(phi) + r2^2 cos^2(psi)
%   Sv : r1^2 cos^2(psi) cos^2(phi) + r2^2 sin^2(psi)
%   Sh : r1^2 sin^2(phi)
%
% so the r1^2 component goes as cos^2(phi) for P and Sv and as sin^2(phi)
% for Sh, while the r2^2 component has no azimuthal dependence at all. The
% weights already in emitP1/kernelCoef are exactly these kernels integrated
% over phi (int cos^2 = int sin^2 = pi, int 1 = 2pi), so adding the azimuth
% refines the direction WITHOUT changing any total: every cross section,
% mean free path and equipartition ratio is untouched.
%   azimuth law   0 : uniform      1 : proportional to cos^2(phi)      2 : to sin^2(phi)
M.emitAzi = [1 0; 1 0; 2 0];              % rows P, Sv, Sh; cols r1^2, r2^2

pg = linspace(0,2*pi,20001).';
M.aziInv = cell(1,2);
cdfa = {(2*pg + sin(2*pg))/(4*pi), ...    % cos^2(phi)
        (2*pg - sin(2*pg))/(4*pi)};       % sin^2(phi)
for c = 1:2
    cc = cdfa{c};
    keep = [true; diff(cc) > 0];          % both laws are flat at their nodes
    M.aziInv{c} = griddedInterpolant(cc(keep),pg(keep),'linear','nearest');
end

% ---- rates (1/s) consumed by the propagator ----------------------------
vB = [alpha0;beta0;beta0];
M.v       = [alpha0;beta0;beta0;cR];
M.rateBB  = n*vB.*M.sigBB;              % 3x3, row = incident mode
M.rateBtot= sum(M.rateBB,2);            % 3x1
M.rateBR  = n*sigBR.*vB.';              % nz x 3, B -> R at depth z
M.rateRB  = n*cR*M.sigRB;               % 1x3
M.rateRR  = n*cR*M.sigRR;
M.rateRtot= sum(M.rateRB) + M.rateRR;

% ---- book-keeping and reference values ---------------------------------
M.alpha0 = alpha0; M.beta0 = beta0; M.freq = freq; M.omega = w;
M.V = V; M.drho = drho; M.n = n; M.zmax = zmax;
M.dz = dz; M.nz = nz; M.zgrid = zg;
M.cR = cR; M.UR = UR; M.kR = kR; M.kP = kP; M.kS = kS;
M.eta_a = ea; M.eta_b = eb; M.r1 = r1; M.r2 = r2;
M.Inorm = Inorm; M.Ir1 = Ir1; M.Ir2 = Ir2; M.secular = secular;

M.tauR = 1/M.rateRtot;                  % Rayleigh mean free time
M.lR   = cR*M.tauR;                     % Rayleigh mean free path
M.lB   = vB./M.rateBtot;                % body mean free paths [P;Sv;Sh]

% Theoretical equipartition ratios Ebar_R / E_B per unit depth (Eq. 17-19)
M.ratioRB = pi*[alpha0 beta0 beta0].^3/(w*cR*UR);
% and among body waves: E_Sv/E_P = E_Sh/E_P = (alpha0/beta0)^3
M.ratioSvP = (alpha0/beta0)^3;

end
