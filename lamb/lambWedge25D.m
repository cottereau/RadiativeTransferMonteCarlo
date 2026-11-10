function out = lambWedge25D(opts)
%LAMBWEDGE25D  Rayleigh wave obliquely incident on a straight 90 degree free wedge (2.5-D SEM).
%   out = lambWedge25D('psi', 30)          psi: angle between the Rayleigh ray and the edge normal
%
%   The edge is the z axis, so the problem is invariant in z: the fields
%       u_x = U_x(x,y) cos(kz z),  u_y = U_y(x,y) cos(kz z),  u_z = U_z(x,y) sin(kz z)
%   close the isotropic elastic equations exactly (a standing pair +-kz; by z-mirror symmetry
%   both carry the same energy fractions). All strain terms then share the z average 1/2, and
%   the cross-section problem is a 2-D SEM with three dofs per node and a kz-dependent K,
%   assembled here from the GLL gradient operators of sem_pml_operators. psi = 0 is the
%   P-SV problem of lambWedgeSEM.
%   kz = k0 sin(psi0) is fixed; each packet component k has psi(k) = asin(kz/k), so the
%   results are averages over psi(k): the band is narrowed so that the energy spread of psi is
%   about dPsi (relWidth = dPsi/tan(psi)), with x0 and rEnd scaled to the packet length
%   (reported: psiRange, the 10-90 % energy range).
%   Cut-offs: P radiates only if vp sin(psi) < cR (psi < 34 deg), S only if vs sin(psi) < cR
%   (psi < 65.5 deg); beyond, the corner only transmits and reflects the Rayleigh wave.
%   Same geometry and measurement as lambWedgeSEM: solid x < 0, y > 0, top face y = 0, side
%   face x = 0, energy snapshot after the corner: R / T strips (2.5 lambda_R), body split into
%   P and S by rho vp^2 div^2 and rho vs^2 |curl|^2, patterns in the in-plane angle
%   theta = atan2(y, -x). Units vs = 1, lambda_S = 1.
arguments
    opts.psi (1,1) double = 0                  % degrees, at the packet centre
    opts.ppw (1,1) double = 8
    opts.vpvs (1,1) double = 4814/2962
    opts.dPsi (1,1) double = 2                 % target spread of psi (deg)
    opts.rEnd (1,1) double = 16                % minimum in-plane distance of the Rayleigh from the corner at tEnd
    opts.nBin (1,1) double = 18
    opts.Yc (1,1) double = 2.5
end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', '..', 'SEMWavePropagation', 'src'));
vs = 1;  vp = opts.vpvs;  rho = 1;
mu = rho*vs^2;  lam = rho*vp^2 - 2*mu;
xi = fzero(@(x) (2 - x)^2 - 4*sqrt(1 - x/vp^2)*sqrt(1 - x), [0.5 0.99]);
cR = sqrt(xi);  qa = sqrt(1 - xi/vp^2);  sa = sqrt(1 - xi);
k0 = 2*pi/cR;  kz = k0*sind(opts.psi);
opts.relWidth = min(0.25, opts.dPsi*pi/180/tand(opts.psi));
sx = cosd(opts.psi)/(k0*opts.relWidth);                            % packet sigma along x
opts.x0 = -(4 + 3.5*sx);  opts.rEnd = max(opts.rEnd, 3.5*sx + 2.5*cR + 8);
cx = cR*cosd(opts.psi);                                            % in-plane (x) group speed of the Rayleigh
tS = opts.rEnd/cx;  tEnd = (-opts.x0)/cx + tS;
vIn = @(c) c*sqrt(max(0, 1 - (c*sind(opts.psi)/cR)^2));          % in-plane group speed of a body wave with this kz (0: cut off)
Rp = max(vIn(vp), vIn(vs))*tS + 2;
p = 4;  h = p/opts.ppw;
L = ceil(max([Rp, -opts.x0 + 3.5*sx + 2, opts.rEnd + 3.5*sx + 4]));   % the whole packet stays off the far faces
mesh = sem_mesh([-L 0; 0 L], round([L L]/h), p);
op = sem_pml_operators(mesh);
nn = size(mesh.nodes, 1);  nq = numel(op.weights);
Q = sparse(1:nq, op.globalNodes, 1, nq, nn);  Z = sparse(nq, nn);  W = spdiags(op.weights, 0, nq, nq);
B1 = op.B{1};  B2 = op.B{2};
Exx = [B1 Z Z];  Eyy = [Z B2 Z];  Ezz = [Z Z kz*Q];
Gxy = [B2 B1 Z];  Gxz = [-kz*Q Z B1];  Gyz = [Z -kz*Q B2];
Dv = Exx + Eyy + Ezz;
K = lam*(Dv'*W*Dv) + 2*mu*(Exx'*W*Exx + Eyy'*W*Eyy + Ezz'*W*Ezz) + mu*(Gxy'*W*Gxy + Gxz'*W*Gxz + Gyz'*W*Gyz);
K = (K + K')/2;
m = rho*full(Q'*op.weights);  Mv = [m; m; m];
clear Exx Eyy Ezz Gxy Gxz Gyz
dtC = sem_critical_dt(K, spdiags(Mv, 0, 3*nn, 3*nn));
fprintf('wedge 2.5-D: psi %g deg (kz/kR %.3f), ppw %g, %d nodes, P radiates %s\n', opts.psi, sind(opts.psi), ...
        opts.ppw, nn, string(vp*sind(opts.psi) < cR));

dt = 0.8*dtC;  nSteps = ceil(tEnd/dt);  dt = tEnd/nSteps;
X = mesh.nodes(:, 1);  Y = mesh.nodes(:, 2);
kk = linspace(max([k0*(1 - 5*opts.relWidth), 0.05*k0, 1.001*kz]), k0*(1 + 5*opts.relWidth), 600);  dk = kk(2) - kk(1);
Ak = exp(-(kk - k0).^2/(2*(opts.relWidth*k0)^2));
u = packet(0);  uPrev = packet(-dt);
Minv = 1./Mv;
E0 = en(u, uPrev);
tic;
for n = 1:nSteps
    uNext = 2*u - uPrev - dt^2*Minv.*(K*u);
    uPrev = u;  u = uNext;
end
uNext = 2*u - uPrev - dt^2*Minv.*(K*u);
fprintf('  %d steps in %.0f s\n', nSteps, toc);
v = (uNext - uPrev)/(2*dt);
E1 = 0.5*sum(Mv.*v.^2) + 0.5*u'*K*u;

% qp densities (all z averages are 1/2: common factor, dropped as in M and K)
U = reshape(u, nn, 3);  V = reshape(v, nn, 3);
ux = U(:, 1);  uy = U(:, 2);  uz = U(:, 3);  qx = Q*ux;  qy = Q*uy;  qz = Q*uz;
exx = B1*ux;  eyy = B2*uy;  ezz = kz*qz;  gxy = B2*ux + B1*uy;  gxz = -kz*qx + B1*uz;  gyz = -kz*qy + B2*uz;
dv = exx + eyy + ezz;
c1 = B2*uz + kz*qy;  c2 = -kz*qx - B1*uz;  c3 = B1*uy - B2*ux;
Vq = Q*V;  w = op.weights;
e = 0.5*rho*sum(Vq.^2, 2) + 0.5*(lam*dv.^2 + 2*mu*(exx.^2 + eyy.^2 + ezz.^2) + mu*(gxy.^2 + gxz.^2 + gyz.^2));
eP = rho*vp^2*dv.^2;  eS = rho*vs^2*(c1.^2 + c2.^2 + c3.^2);
xq = mesh.nodes(op.globalNodes, 1);  yq = mesh.nodes(op.globalNodes, 2);
Yc = opts.Yc*cR;
inR = yq < Yc & xq < -Yc;  inT = xq > -Yc & yq > Yc;  body = ~inR & ~inT;
psiK = asind(kz./kk);  Ek = Ak.^2.*kk.^2;  cE = cumsum(Ek)/sum(Ek);
out = struct('psi', opts.psi, 'ppw', opts.ppw, 'psiRange', [interp1(cE, psiK, 0.9), interp1(cE, psiK, 0.1)], ...
             'conservation', E1/E0 - 1, 'qpVsMK', sum(w.*e)/E1 - 1);
out.R = sum(w(inR).*e(inR))/E0;  out.T = sum(w(inT).*e(inT))/E0;  out.body = sum(w(body).*e(body))/E0;
out.P = sum(w(body).*eP(body))/E0;  out.S = sum(w(body).*eS(body))/E0;
out.bodySplitCheck = (out.P + out.S)/out.body - 1;  out.sum = out.R + out.T + out.body;
th = atan2(yq, -xq);  ed = linspace(0, pi/2, opts.nBin + 1);  out.thetaEdges = ed;
b = discretize(th(body), ed);  ok = ~isnan(b);  wb = w(body);  bP = eP(body);  bS = eS(body);
out.patP = accumarray(b(ok), wb(ok).*bP(ok), [opts.nBin 1]).'/E0;
out.patS = accumarray(b(ok), wb(ok).*bS(ok), [opts.nBin 1]).'/E0;
fprintf('  psi %.1f-%.1f: E conservation %.1e, qp vs M/K %.1e\n', out.psiRange, out.conservation, out.qpVsMK);
fprintf('  R %.4f  T %.4f  body %.4f (P %.4f, S %.4f, split check %.1e)  sum %.4f\n', ...
        out.R, out.T, out.body, out.P, out.S, out.bodySplitCheck, out.sum);

    function uu = packet(t)                                        % right-going Rayleigh, wavevector (kx, kz)
        Ux = zeros(nn, 1);  Uy = Ux;  Uz = Ux;
        for c = 1:20:numel(kk)
            kc = kk(c:min(c + 19, end));  kx = sqrt(kc.^2 - kz^2);
            ph = exp(1i*(kx.*(X - opts.x0) - cR*kc*t));
            eq = exp(-qa*Y*kc);  es = exp(-sa*Y*kc);
            ul = ph.*(1i*(eq - 2*qa*sa/(sa^2 + 1)*es));             % along the ray, per unit k
            uv = ph.*(-qa*(eq - 2/(sa^2 + 1)*es));
            a = Ak(c:min(c + 19, end))*dk;
            Ux = Ux + real(ul*(a.*kx./kc).');
            Uy = Uy + real(uv*a.');
            Uz = Uz + real(1i*ul*(a*kz./kc).');                     % sin(kz z) partner of the +-kz pair
        end
        uu = [Ux; Uy; Uz];
    end

    function E = en(a, aPrev)
        aNext = 2*a - aPrev - dt^2*Minv.*(K*a);
        vv = (aNext - aPrev)/(2*dt);
        E = 0.5*sum(Mv.*vv.^2) + 0.5*a'*K*a;
    end
end
