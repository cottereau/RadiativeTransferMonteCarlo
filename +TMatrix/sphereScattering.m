function out = sphereScattering(a,m1,m2,omega,opts)
%% sphereScattering - elastic scattering by a single elastic sphere
%
%   out = TMatrix.sphereScattering(a, host, inclusion, omega)
%
% host and inclusion: structs with fields lam, mu, rho.
%
% P incidence is axisymmetric (m=0) and gives one 4x4 system per order n.
% S incidence lives in m=1 and splits into a 4x4 spheroidal system (P + SV)
% plus a 2x2 toroidal one (SH), decoupled by orthogonality of the vector harmonics.
%
% The incident plane-wave coefficients are NOT written from memory: they are
% extracted by numerical projection of the incident field onto the basis. The
% same holds for the far-field amplitudes, obtained by evaluating the field at a large radius.
%
% Output:
%   out.th                   scattering angles (rad)
%   out.dsig_PP, dsig_PS     differential cross-section for P incidence (m^2/sr)
%   out.dsig_SP, dsig_SS     same for S incidence, averaged over azimuth
%                            and over incident polarization
%   out.sig_PP ... sig_SS    total cross-sections (m^2)
%   out.nmax

if nargin < 5, opts = struct; end
if ~isfield(opts,'nth'),  opts.nth  = 361;  end
if ~isfield(opts,'rfar'), opts.rfar = 2000; end
if ~isfield(opts,'nmax'), opts.nmax = [];   end
if ~isfield(opts,'nquad'),opts.nquad= 2000; end
if ~isfield(opts,'inner'),opts.inner= 0;   end   % B21: nr > 0 -> internal field on nr radii (out.Win)

vp1 = sqrt((m1.lam+2*m1.mu)/m1.rho);  vs1 = sqrt(m1.mu/m1.rho);
vp2 = sqrt((m2.lam+2*m2.mu)/m2.rho);  vs2 = sqrt(m2.mu/m2.rho);
kp1 = omega/vp1; ks1 = omega/vs1;
kp2 = omega/vp2; ks2 = omega/vs2;

if isempty(opts.nmax)
    xx = max(ks1,ks2)*a;
    opts.nmax = ceil(xx + 4*xx^(1/3) + 6);
end
N  = opts.nmax;
rf = opts.rfar*a;

th = linspace(1e-6,pi-1e-6,opts.nth);
u  = cos(th);

% ---- projection of the incident plane wave onto the basis (numerical) ----
inc = projectIncident(N,ks1,kp1,a,m1,opts.nquad);

% internal field (B21, opts.inner = nr): components in (r, th), r at the middle of nr shells
ri = ((1:opts.inner)' - 0.5)/max(opts.inner,1)*a;
iP = zeros(numel(ri),numel(th),2); iS = zeros(numel(ri),numel(th),3);
% accumulated far fields
urP = zeros(size(th)); uthP = zeros(size(th));                    % incid. P
urS = zeros(size(th)); uthS = zeros(size(th)); uphS = zeros(size(th)); % incid. S

for n = 0:N
    % radial factors at r=a
    Ph_a  = radial('h',n,0,kp1,a,m1,'P');    % m=0
    Sh_a  = radial('h',n,0,ks1,a,m1,'SV');
    Pj2_a = radial('j',n,0,kp2,a,m2,'P');
    Sj2_a = radial('j',n,0,ks2,a,m2,'SV');
    Pj1_a = radial('j',n,0,kp1,a,m1,'P');

    % ---------- P incidence (m=0) ----------
    lin = [1 2 4 5];                        % ur, uth, srr, srt
    if n == 0
        M  = [Ph_a([1 4]), -Pj2_a([1 4])];
        xx = M \ (-inc.aP(n+1)*Pj1_a([1 4]));
        A0 = xx(1); B0 = 0;  C0 = xx(2); D0 = 0;
    else
        M  = [Ph_a(lin), Sh_a(lin), -Pj2_a(lin), -Sj2_a(lin)];
        M  = balanceColumns(M);
        xx = M.M \ (-inc.aP(n+1)*Pj1_a(lin));
        xx = xx./M.s;
        A0 = xx(1); B0 = xx(2);  C0 = xx(3); D0 = xx(4);
    end

    % ---------- S incidence (m=1) ----------
    A1 = 0; B1 = 0; E1 = 0;
    if n >= 1
        Ph1  = radial('h',n,1,kp1,a,m1,'P');
        Sh1  = radial('h',n,1,ks1,a,m1,'SV');
        Pj21 = radial('j',n,1,kp2,a,m2,'P');
        Sj21 = radial('j',n,1,ks2,a,m2,'SV');
        Sj11 = radial('j',n,1,ks1,a,m1,'SV');

        M  = [Ph1(lin), Sh1(lin), -Pj21(lin), -Sj21(lin)];
        M  = balanceColumns(M);
        xx = M.M \ (-inc.cS(n+1)*Sj11(lin));
        xx = xx./M.s;
        A1 = xx(1); B1 = xx(2);  C1 = xx(3); D1 = xx(4);

        % toroidal: 2x2 in (uth, srt)
        Th1  = radial('h',n,1,ks1,a,m1,'T');
        Tj21 = radial('j',n,1,ks2,a,m2,'T');
        Tj11 = radial('j',n,1,ks1,a,m1,'T');
        lt = [2 5];
        Mt = balanceColumns([Th1(lt), -Tj21(lt)]);
        xt = Mt.M \ (-inc.dS(n+1)*Tj11(lt));
        xt = xt./Mt.s;
        E1 = xt(1);  F1 = xt(2);
    end

    % ---------- far field ----------
    [P0,P1,~,~] = TMatrix.sphereFun('P',n,u);
    dP0 = -sin(th).*P1;                        % dP_n/dth
    if n >= 1
        [Q0,Q1,~,~] = TMatrix.legAssoc(n,1,u);
        dQ0 = -sin(th).*Q1;                    % dP_n^1/dth
        aQ  = Q0./sin(th);                     % (m/sin th) P_n^1, m=1
    end

    gP  = farFat('P', n,0,kp1,rf,m1,1);        % ur of the P field
    gSV = farFat('SV',n,0,ks1,rf,m1,2);        % uth of the SV field

    for ir = 1:numel(ri)                       % B21: internal field
        Pr = radial('j',n,0,kp2,ri(ir),m2,'P');  Sr = radial('j',n,0,ks2,ri(ir),m2,'SV');
        iP(ir,:,1) = iP(ir,:,1) + (C0*Pr(1) + D0*Sr(1))*P0;
        iP(ir,:,2) = iP(ir,:,2) + (C0*Pr(2) + D0*Sr(2))*dP0;
        if n >= 1
            Pr = radial('j',n,1,kp2,ri(ir),m2,'P');  Sr = radial('j',n,1,ks2,ri(ir),m2,'SV');  Tr = radial('j',n,1,ks2,ri(ir),m2,'T');
            iS(ir,:,1) = iS(ir,:,1) + (C1*Pr(1) + D1*Sr(1))*Q0;
            iS(ir,:,2) = iS(ir,:,2) + (C1*Pr(2) + D1*Sr(2))*dQ0 + F1*Tr(2)*aQ;
            iS(ir,:,3) = iS(ir,:,3) + (C1*Pr(3) + D1*Sr(3))*aQ  + F1*Tr(3)*dQ0;
        end
    end
    urP  = urP  + A0*gP *P0;
    uthP = uthP + B0*gSV*dP0;

    if n >= 1
        gP1  = farFat('P', n,1,kp1,rf,m1,1);   % ur
        gSt  = farFat('SV',n,1,ks1,rf,m1,2);   % uth
        gSp  = farFat('SV',n,1,ks1,rf,m1,3);   % uph
        gTt  = farFat('T', n,1,ks1,rf,m1,2);
        gTp  = farFat('T', n,1,ks1,rf,m1,3);

        urS  = urS  + A1*gP1*Q0;
        uthS = uthS + B1*gSt*dQ0 + E1*gTt*aQ;
        uphS = uphS + B1*gSp*aQ  + E1*gTp*dQ0;
    end
end

% ---- differential cross-sections ----
% P incidence: no azimuthal dependence
dsig_PP = abs(urP).^2;
dsig_PS = abs(uthP).^2*(vs1/vp1);

% S incidence: azimuthal average (cos^2 and sin^2 -> 1/2), which also averages
% over the incident polarization
dsig_SP = (abs(urS).^2/2)*(vp1/vs1);
dsig_SS = (abs(uthS).^2 + abs(uphS).^2)/2;

w = 2*pi*sin(th);
out = struct('th',th, ...
    'dsig_PP',dsig_PP,'dsig_PS',dsig_PS,'dsig_SP',dsig_SP,'dsig_SS',dsig_SS, ...
    'sig_PP',trapz(th,w.*dsig_PP),'sig_PS',trapz(th,w.*dsig_PS), ...
    'sig_SP',trapz(th,w.*dsig_SP),'sig_SS',trapz(th,w.*dsig_SS), ...
    'nmax',N,'kp',kp1,'ks',ks1,'a',a);

% COMPLEX amplitudes: they keep the phase, needed for the optical theorem.
% As theta -> 0 the scattered field becomes parallel to the incident polarization,
% so f(0) is the dominant component evaluated at the first grid node.
out.f_PP = urP;  out.f_PS = uthP;
out.f_SS = uthS; out.f_SP = urS;
out.f_SSph = uphS;             % phi component of the SS field (B21 Wigner delay)
out.f_PP_fwd = urP(1);          % e_r -> e_z at theta=0
out.f_SS_fwd = uthS(1);         % e_th -> e_x at theta=0, phi=0
out.sig_P = out.sig_PP + out.sig_PS;
if opts.inner > 0      % B21: int |u|^2 dV of the internal field, P and S incidence (S averaged in azimuth: cos^2, sin^2 -> 1/2)
    vol = @(F) 2*pi*sum(ri.^2.*trapz(th, sin(th).*F, 2))*(a/opts.inner);
    out.U2in = [vol(sum(abs(iP).^2,3)), vol(sum(abs(iS).^2,3)/2)];
end
out.sig_S = out.sig_SS + out.sig_SP;
end

% =======================================================================
function Q = radial(kind,n,m,k,r,med,type)
% six radial factors [ur uth uph srr srt srp], already divided by the matching
% angular function. The factor is independent of th and phi (validated in
% testSphereBCm), so any point will do.
th = 0.7; ph = 0.4;
[Z0,Z1,Z2,Z3] = TMatrix.sphereFun(kind,n,k*r);
Zr = [Z0, k*Z1, k^2*Z2, k^3*Z3];
[L0,L1,~,~]   = TMatrix.legAssoc(n,m,cos(th));
[~,~,L2,L3]   = TMatrix.legAssoc(n,m,cos(th));
bc = TMatrix.sphereBCkernel(r,th,ph,m,med.lam,med.mu,Zr(1),Zr(2),Zr(3),Zr(4),L0,L1,L2,L3);

C = cos(m*ph); S = sin(m*ph);
dP = -sin(th)*L1;
ax = m*L0/sin(th);

switch type
    case 'P',  b = bc(1:6);   ang = [L0*C; dP*C; ax*S; L0*C; dP*C; ax*S];
    case 'SV', b = bc(7:12);  ang = [L0*C; dP*C; ax*S; L0*C; dP*C; ax*S];
    case 'T',  b = bc(13:18); ang = [1;    ax*C; dP*S; 1;    ax*C; dP*S];
end
Q = b./ang;
Q(~isfinite(Q)) = 0;      % components that vanish by symmetry (m=0)
end

function g = farFat(type,n,m,k,rf,med,comp)
% value of component 'comp' of the far field per unit coefficient, with the
% outgoing-wave factor and the 1/r removed
Q = radial('h',n,m,k,rf,med,type);
g = Q(comp)*rf*exp(-1i*k*rf);
end

function E = balanceColumns(M)
% scale the columns to unit norm before solving: otherwise the system has
% RCOND ~1e-21 in the ka->0 limit, because the four basic solutions have
% very different scales
s = vecnorm(M,2,1);
s(s==0) = 1;
E.M = M./s;
E.s = s(:);
end

% =======================================================================
function inc = projectIncident(N,ks,kp,a,med,nq)
% Incident plane-wave coefficients in the partial-wave basis, obtained by
% numerical projection -- no expansion constant written from memory.
%
%   P wave:  u = e_z exp(i kp z)   -> spheroidal m=0 only
%   S wave:  u = e_x exp(i ks z)   -> spheroidal m=1 (cos ph) + toroidal m=1
%
% Project u_r (which only receives spheroidal contributions) to get aP and cS,
% then u_ph to get dS by least squares.

thq = linspace(1e-6,pi-1e-6,nq);
wq  = sin(thq);
uq  = cos(thq);

inc.aP = zeros(N+1,1);
inc.cS = zeros(N+1,1);
inc.dS = zeros(N+1,1);

% --- radial component of the incident fields ---
urP_inc = cos(thq).*exp(1i*kp*a*cos(thq));           % e_z . r_hat
urS_inc = sin(thq).*exp(1i*ks*a*cos(thq));           % coefficient of cos(ph) in e_x . r_hat
uphS_inc= -exp(1i*ks*a*cos(thq));                    % coefficient of sin(ph) in e_x . ph_hat

for n = 0:N
    [P0,~,~,~] = TMatrix.sphereFun('P',n,uq);
    Rp = radial('j',n,0,kp,a,med,'P');
    nrm = trapz(thq,wq.*P0.^2);
    if Rp(1) ~= 0 && nrm > 0
        inc.aP(n+1) = trapz(thq,wq.*urP_inc.*P0)/(nrm*Rp(1));
    end

    if n >= 1
        [Q0,~,~,~] = TMatrix.legAssoc(n,1,uq);
        Rs = radial('j',n,1,ks,a,med,'SV');
        nrmQ = trapz(thq,wq.*Q0.^2);
        if Rs(1) ~= 0 && nrmQ > 0
            inc.cS(n+1) = trapz(thq,wq.*urS_inc.*Q0)/(nrmQ*Rs(1));
        end
    end
end

% --- toroidal part by least squares on u_ph ---
% u_ph(sin ph) = sum_n [ cS_n * Rs(3) * (P_n^1/sin th) + dS_n * Rt(3) * dP_n^1/dth ]
A = zeros(nq,N); rhs = uphS_inc(:);
for n = 1:N
    [Q0,Q1,~,~] = TMatrix.legAssoc(n,1,uq);
    dQ = -sin(thq).*Q1;
    aQ = Q0./sin(thq);
    Rs = radial('j',n,1,ks,a,med,'SV');
    Rt = radial('j',n,1,ks,a,med,'T');
    rhs = rhs - inc.cS(n+1)*Rs(3)*aQ(:);
    A(:,n) = Rt(3)*dQ(:);
end
W = sqrt(wq(:));
inc.dS(2:N+1) = (W.*A) \ (W.*rhs);
end
