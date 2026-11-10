function P = propagateParticleRayleigh(M,P,dt,fsbc)
%PROPAGATEPARTICLERAYLEIGH  One sub-step of the coupled body/Rayleigh walk.
%
%   P = propagateParticleRayleigh(M,P,dt,fsbc) advances the particle
%   ensemble P by one time sub-step dt under the model M built by
%   RayleighModel.  fsbc selects the free surface boundary condition,
%   'modeconv' (traction free, P<->Sv conversion) or 'mirror' (specular).
%   Xu et al. (2021, Table 1) show both give the same equipartition ratios.
%
%   P has fields
%     x    : N x 3 positions [x y z], z positive downwards, surface at z = 0
%     dir  : N x 3 unit propagation directions (dir(:,3) = 0 for Rayleigh)
%     mode : N x 1 in {1,2,3,4} = {P, Sv, Sh, Rayleigh}
%
%   The depth-dependent body-to-Rayleigh coupling is integrated exactly
%   along the step (via the tabulated vertical optical depth M.FBR) rather
%   than sampled at the step ends, so dt is set by the position resolution
%   wanted, not by the 0.5 km scale over which sigma_BR(z) decays.
%
%   The scatterers fill the whole halfspace: M.zmax only truncates the
%   Rayleigh coupling tables, where sigma_BR has already decayed. Confining
%   the scatterers to 0 < z < zmax instead would let every body wave escape
%   downwards (l_P = 147 km against a 10 km layer) and no equipartition
%   would ever set in.
%
%   See also RAYLEIGHMODEL, VALIDATIONXU2021.
%
%   References:
%     - Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Monte Carlo simulations of coupled
%       body- and Rayleigh-wave multiple scattering in elastic media, Geophys. J. Int. 228,
%       1213-1236. doi:10.1093/gji/ggab370

N  = numel(P.mode);
z0 = P.x(:,3);
m0  = P.mode;                           % state before the surface interaction
mu0 = P.dir(:,3);

% ---- straight-line advance ---------------------------------------------
P.x = P.x + (M.v(P.mode)*dt).*P.dir;

% ---- free surface at z = 0 ---------------------------------------------
hit = (P.mode < 4) & (P.x(:,3) < 0);
if any(hit)
    P.x(hit,3) = -P.x(hit,3);
    if strcmpi(fsbc,'modeconv')
        [P.dir(hit,:),P.mode(hit)] = ...
            freeSurface(P.dir(hit,:),P.mode(hit),M.alpha0,M.beta0);
    else
        P.dir(hit,3) = -P.dir(hit,3);
    end
end

m   = P.mode;                           % state after the surface interaction
mu  = P.dir(:,3);
z1  = P.x(:,3);
B   = m < 4;
col = min(m,3);

% ---- optical depth accumulated over the step ---------------------------
% body -> Rayleigh. The kernel (B.3-B.5) weights the two depth profiles by
% the incidence angle, sigma_BR = a(mu)*r1^2(z) + b(mu)*r2^2(z), and each
% profile integrates exactly along the straight step as |G(z1)-G(z0)|/|mu|
% (two legs z0 -> 0 -> z1 when the surface is hit, G(0) = 0).
% A step that hits the surface is split into its two legs, because the
% traction-free surface converts most incident P into Sv and the legs then
% have different kernels, speeds and mu -- right where sigma_BR peaks.
[a0,b0] = kernelCoef(M,m0,mu0);
G10 = tabLookup(M.G1,M,z0);  G20 = tabLookup(M.G2,M,z0);
G11 = tabLookup(M.G1,M,z1);  G21 = tabLookup(M.G2,M,z1);
amu0 = abs(mu0);

tauBR = (a0.*abs(G11-G10) + b0.*abs(G21-G20))./max(amu0,1e-12);
if any(hit)
    [a1,b1] = kernelCoef(M,m(hit),mu(hit));
    amu1 = abs(mu(hit));
    tauBR(hit) = (a0(hit).*G10(hit) + b0(hit).*G20(hit))./max(amu0(hit),1e-12) ...
               + (a1     .*G11(hit) + b1     .*G21(hit))./max(amu1,1e-12);
end
flat = B & (amu0 < 1e-5);               % near horizontal: use the local rate
if any(flat)
    zl = z0(flat);
    tauBR(flat) = (a0(flat).*M.r1(zl).^2 + b0(flat).*M.r2(zl).^2) ...
                  .*(M.n*M.v(m0(flat))*dt);
end
tauBR(~B) = 0;

tau = zeros(N,1);
tau(B)  = M.rateBtot(col(B))*dt + tauBR(B);
tau(~B) = M.rateRtot*dt;

% ---- scattering events -------------------------------------------------
k = find(rand(N,1) < -expm1(-tau));
if isempty(k), return; end
nk = numel(k);
Bk = B(k);

W = zeros(nk,4);
if any(Bk)
    kk = k(Bk);
    W(Bk,1:3) = M.rateBB(col(kk),:)*dt;
    W(Bk,4)   = tauBR(kk);
end
if any(~Bk)
    W(~Bk,1:3) = repmat(M.rateRB*dt,nnz(~Bk),1);
    W(~Bk,4)   = M.rateRR*dt;
end
C   = cumsum(W,2)./sum(W,2);
tgt = sum(rand(nk,1) > C(:,1:3),2) + 1;

fromR = ~Bk;
toR   = tgt == 4;

% Rayleigh -> body: reciprocity makes the outgoing (depth, colatitude) pair
% follow the same kernel, a two-component mixture of z ~ r1^2 and z ~ r2^2,
% each with its own law for mu = cos(psi1) (Appendices D and E).
% The azimuth follows B.3-B.5 too: the r1^2 component goes as cos^2(phi)
% (P, Sv) or sin^2(phi) (Sh), the r2^2 component is uniform. phi is measured
% from the direction the Rayleigh wave was travelling in. Sampling it changes
% no total -- the weights in emitP1 are these kernels already integrated over
% phi -- it only stops the outgoing body wave from being azimuthally flat.
j = fromR & ~toR;
if any(j)
    kj = k(j);  tj = tgt(j);  nj = numel(kj);
    cp = 1 + (rand(nj,1) > M.emitP1(tj));      % 1: z ~ r1^2,  2: z ~ r2^2
    u  = rand(nj,1);
    zz = zeros(nj,1);
    for c = 1:2
        sc = cp == c;
        if any(sc), zz(sc) = M.zInv{c}(u(sc)); end
    end
    law = M.emitLaw(sub2ind(size(M.emitLaw),tj,cp));
    ur  = rand(nj,1);  mm = zeros(nj,1);
    sc = law==0;  mm(sc) = 2*ur(sc)-1;                            % uniform
    sc = law==1;  mm(sc) = M.muInv(ur(sc));                       % ~ 1-mu^2
    sc = law==2;  mm(sc) = sign(2*ur(sc)-1).*abs(2*ur(sc)-1).^(1/3); % ~ mu^2
    % azimuth relative to the incoming Rayleigh direction
    aziLaw = M.emitAzi(sub2ind(size(M.emitAzi),tj,cp));
    ur2 = rand(nj,1);  dphi = 2*pi*ur2;
    for c = 1:2
        sc = aziLaw == c;
        if any(sc), dphi(sc) = M.aziInv{c}(ur2(sc)); end
    end
    psiR = atan2(P.dir(kj,2), P.dir(kj,1));   % Rayleigh travels horizontally
    ph = psiR + dphi;
    st = sqrt(1-mm.^2);
    P.x(kj,3)   = zz;
    P.dir(kj,:) = [st.*cos(ph) st.*sin(ph) mm];
end

% body -> Rayleigh: the particle collapses onto the free surface
P.x(k(~fromR & toR),3) = 0;

% ---- outgoing directions -----------------------------------------------
% ponytail: body-to-body scattering is isotropic, using the angle-averaged
% totals B.36-B.39. Detailed balance fixes the equipartition ratios from
% the totals alone; add the phase functions B.29-B.33 if transport mean
% free paths or angular distributions are ever needed.
kb = k(~toR & ~fromR);
if ~isempty(kb)
    c  = 2*rand(numel(kb),1)-1;
    ph = 2*pi*rand(numel(kb),1);
    s  = sqrt(1-c.^2);
    P.dir(kb,:) = [s.*cos(ph) s.*sin(ph) c];
end
kr = k(toR);
if ~isempty(kr)
    % R -> R keeps a uniform azimuth: B.24 gives the total only, and neither
    % Xu's appendix nor the modal derivation supplies an R->R phase function.
    ph = 2*pi*rand(numel(kr),1);
    isB2R = ~fromR(toR);                   % body -> Rayleigh, B.3-B.5
    if any(isB2R)
        ph(isB2R) = bodyToRayleighAzimuth(M, P, kr(isB2R), m);
    end
    P.dir(kr,:) = [cos(ph) sin(ph) zeros(numel(kr),1)];
end
P.mode(k) = tgt;

end

% -------------------------------------------------------------------------
function ph = bodyToRayleighAzimuth(M, P, kr, m)
% Azimuth of the Rayleigh wave a body wave converts into, from B.3-B.5.
%
%   P  : r1^2 sin^2(psi) cos^2(phi) + r2^2 cos^2(psi)
%   Sv : r1^2 cos^2(psi) cos^2(phi) + r2^2 sin^2(psi)
%   Sh : r1^2 sin^2(phi)
%
% kernelCoef already returns a(mu) and b(mu), which ARE these two kernels
% integrated over phi, so a*r1^2 and b*r2^2 are the component weights and no
% new normalisation is needed. phi is measured from the horizontal direction
% of the incident ray; for a ray arriving vertically that direction is
% degenerate, and the polarization supplies it instead (for P the r1^2 term
% vanishes there anyway, so the azimuth is uniform and no axis is needed).
n = numel(kr);
mo = m(kr);
mu = P.dir(kr,3);
[a,b] = kernelCoef(M, mo, mu);
r1z = M.r1(max(P.x(kr,3),0));
r2z = M.r2(max(P.x(kr,3),0));

w1 = a.*r1z.^2;                       % r1^2 component (phi-dependent)
w2 = b.*r2z.^2;                       % r2^2 component (uniform in phi)
tot = w1 + w2;
comp1 = rand(n,1) .* max(tot,realmin) < w1;

% reference axis: horizontal projection of the ray, else of the polarization
ax = P.dir(kr,1:2);
bad = vecnorm(ax,2,2) < 1e-9;
if any(bad)
    alt = P.perp(kr,1:2);
    ax(bad,:) = alt(bad,:);
end
na = vecnorm(ax,2,2);
degenerate = na < 1e-9;               % measure zero: fall back to uniform
na(degenerate) = 1;
psi0 = atan2(ax(:,2)./na, ax(:,1)./na);
psi0(degenerate) = 2*pi*rand(nnz(degenerate),1);

ur = rand(n,1);
dphi = 2*pi*ur;                       % r2^2 component, and the degenerate rows
lawC = M.emitAzi(sub2ind(size(M.emitAzi), mo, ones(n,1)));   % r1^2 law per mode
for c = 1:2
    sc = comp1 & lawC == c & ~degenerate;
    if any(sc), dphi(sc) = M.aziInv{c}(ur(sc)); end
end
ph = psi0 + dphi;
end

% -------------------------------------------------------------------------
function [a,b] = kernelCoef(M,m,mu)
% sigma_BR = a(mu)*r1^2(z) + b(mu)*r2^2(z), from B.3-B.5
mu2 = mu.^2;
a = zeros(numel(m),1);  b = zeros(numel(m),1);
i1 = m==1;  i2 = m==2;  i3 = m==3;
a(i1) = M.cP*(1-mu2(i1));   b(i1) = 2*M.cP*mu2(i1);
a(i2) = M.cS*mu2(i2);       b(i2) = 2*M.cS*(1-mu2(i2));
a(i3) = M.cS;
end

% -------------------------------------------------------------------------
function y = tabLookup(tab,M,z)
% linear interpolation on the uniform depth grid, saturating at zmax
t = min(max(z,0),M.zmax)/M.dz;
i = min(floor(t)+1,M.nz-1);
w = t-(i-1);
y = tab(i).*(1-w) + tab(i+1).*w;
end

% -------------------------------------------------------------------------
function [d,mo] = freeSurface(d,mo,a,b)
% Traction-free surface. Sh reflects specularly; for P and Sv the surface
% reflects all the energy, so the converted share is exactly 1 - R^2 with R
% the unconverted displacement reflection coefficient (Aki & Richards 5.26).
si = sqrt(max(0,1-d(:,3).^2));          % sine of the incidence angle
hu = d(:,1:2)./max(si,realmin);         % horizontal unit vector
d(:,3) = -d(:,3);                       % specular, now pointing downwards

for inc = 1:2                           % 1 = P incident, 2 = Sv incident
    q = find(mo == inc);
    if isempty(q), continue; end
    if inc == 1
        p = si(q)/a;  vin = a;  vout = b;  mout = 2;
    else
        p = si(q)/b;  vin = b;  vout = a;  mout = 1;
    end
    qa = sqrt(max(0,1/a^2-p.^2));
    qb = sqrt(max(0,1/b^2-p.^2));
    A  = (1/b^2-2*p.^2).^2;
    Bc = 4*p.^2.*qa.*qb;
    R  = (A-Bc)./(A+Bc);                % beyond the Sv critical angle Bc = 0
    cv = rand(numel(q),1) > R.^2;
    if ~any(cv), continue; end
    qc = q(cv);
    sn = min(1,si(qc)*vout/vin);        % Snell's law
    d(qc,:) = [hu(qc,:).*sn , sqrt(max(0,1-sn.^2))];
    mo(qc)  = mout;
end
end
