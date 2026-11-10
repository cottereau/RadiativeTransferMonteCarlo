function C = interfaceLameChannels(K, zlim, w, vp, vs, rho, rhoRef, cI, UI, Inorm, V, dlam, dmu)
%INTERFACELAMECHANNELS  delta-lambda and delta-mu channels of an interface mode.
%
%   C = interfaceLameChannels(...) returns, for heterogeneity in the Lame
%   parameters inside ONE medium adjoining the interface, both directions of
%   coupling to that medium's bulk waves: the depth-integrated
%   interface -> bulk cross sections, and the depth-resolved
%   bulk -> interface coefficients.
%
%   These generalise Xu et al. (2021) B.1/B.15 (delta-lambda) and
%   B.12-B.14 / B.20-B.22 (delta-mu) from a free surface to an interface,
%   by the same rule already validated for delta-rho: the local material
%   constants become medium j's, and the whole channel carries the weight
%
%       rho_j / rho_ref
%
%   that the shared modal normalisation Inorm = Inorm_1 + (rho_2/rho_1)Inorm_2
%   demands. That weight is not free -- it cancels against itself in
%   detailed balance, leaving the medium-independent ratio pi v^3/(w c U)
%   which testInterfaceWavesAndWedge asserts for every channel.
%
%   A fluid has mu = 0, so pass dmu = 0 there: only the delta-lambda channel
%   exists, and only into the acoustic wave.
%
%   C.sig_lam_P, C.sig_mu_P, C.sig_mu_Sv, C.sig_mu_Sh   interface -> bulk (m^2)
%   C.c_lam_P,   C.c_mu_P,   C.c_mu_Sv,   C.c_mu_Sh     bulk -> interface (m^4)
%   C.IK_lam,    C.IK_mu_P,  C.IK_mu_Sv,  C.IK_mu_Sh    depth integrals
%
%   See also INTERFACELAMEKERNELS, SCHOLTEMODEL, STONELEYMODEL.
%
%   References:
%     - Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Monte Carlo simulations of coupled
%       body- and Rayleigh-wave multiple scattering in elastic media, Geophys. J. Int. 228,
%       1213-1236. doi:10.1093/gji/ggab370

mu  = rho*vs^2;
lam = rho*vp^2 - 2*mu;
kP  = w/vp;
wj  = rho/rhoRef;
e2lam = V^2*dlam^2;
e2mu  = V^2*dmu^2;

I = @(f) integral(f, zlim(1), zlim(2), 'RelTol', 1e-9, 'AbsTol', 0);
C.IK_lam   = I(K.lam_P);
if vs > 0
    kS = w/vs;
    C.IK_mu_P  = I(K.mu_P);
    C.IK_mu_Sv = I(K.mu_Sv);
    C.IK_mu_Sh = I(K.mu_Sh);
else
    kS = 0;
    C.IK_mu_P = 0;  C.IK_mu_Sv = 0;  C.IK_mu_Sh = 0;
end

% ---- interface -> bulk, already integrated over the scatterer depth ----
C.sig_lam_P  = kP^4*(lam/rho)^2*wj/(4*pi*vp*cI^2*UI*Inorm) * C.IK_lam * e2lam;
if vs > 0
    C.sig_mu_P  = kP^4*(mu/rho)^2*wj/(15*pi*vp*cI^2*UI*Inorm) * C.IK_mu_P  * e2mu;
    C.sig_mu_Sv = w^3*kS*wj/(pi*cI^2*UI*Inorm)                * C.IK_mu_Sv * e2mu;
    C.sig_mu_Sh = w^3*kS*wj/(pi*cI^2*UI*Inorm)                * C.IK_mu_Sh * e2mu;
else
    C.sig_mu_P = 0;  C.sig_mu_Sv = 0;  C.sig_mu_Sh = 0;
end

% ---- bulk -> interface, per unit depth: multiply by the kernel at z_s ---
C.c_lam_P  = kP^3*(lam/rho)^2*wj/(4*cI*UI^3*Inorm) * e2lam;
if vs > 0
    C.c_mu_P  = kP^3*(mu/rho)^2*wj/(15*cI*UI^3*Inorm) * e2mu;
    C.c_mu_Sv = w^3*vs*wj/(cI*UI^3*Inorm)             * e2mu;
    C.c_mu_Sh = w^3*vs*wj/(cI*UI^3*Inorm)             * e2mu;
else
    C.c_mu_P = 0;  C.c_mu_Sv = 0;  C.c_mu_Sh = 0;
end
end
