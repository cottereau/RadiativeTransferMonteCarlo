function out = lambWedgeSEM(opts)
%LAMBWEDGESEM  Rayleigh wave at a 90 degree free wedge, 2-D P-SV SEM.
%   out = lambWedgeSEM('ppw', 8)                 corner run
%   out = lambWedgeSEM('ppw', 8, 'corner', false) reference half-space (strip retention check)
%   out = lambWedgeSEM('incident', 'P', 'theta', 45)  inverse path: P (or S) beam from angle theta
%
%   Units vs = 1, lambda_S = 1 (the wedge has no length scale: T, R, P, S depend only on
%   vp/vs; ppw, the GLL intervals per lambda_S at the packet centre, plays the role of the
%   frequency: 250 / 500 / 750 kHz on one mesh = ppw x 2 / 1 / 2/3).
%   Solid x < 0, y > 0 (y = depth); free faces y = 0 (top) and x = 0 (side), the corner at the
%   origin. The far faces are free too but nothing reaches them before tEnd.
%   Initial condition: an exact right-going Rayleigh packet (Gaussian in k, centre k_R) on the
%   top face at x0, so no body waves are launched. After the corner (Rayleigh at distance
%   rEnd from it) the energy is split at qp level:
%     R: top strip  y < Yc, x < -Yc;   T: side strip  x > -Yc, y > Yc;   body: the rest,
%   Yc = 2.5 lambda_R; body into P and S by the plane-wave densities rho vp^2 (div u)^2 and
%   rho vs^2 |curl u|^2, checked against the body energy density. Angular patterns of P and S
%   in theta = atan2(y, -x) (0 = along the top face, 90 deg = down the side face).
%   Inverse path: a Gaussian beam (width w, energy profile exp(-m^2/w^2)) of plane P or S waves
%   coming from angle theta; the Rayleigh energy W on each face (strips beyond rEnd/2 from the
%   corner) over the fluence F at the beam waist is sigma = W/F; the spreading of the beam to the
%   corner (z_R = k w^2) is folded into <lambda_c> below.
%   Detailed balance (equipartition: 2-D body states k dk dtheta/(2 pi)^2, Rayleigh dk/(2 pi)
%   per face and direction) with reciprocity gives sigma_{c,theta -> top} = <lambda_c> dE_c/dtheta
%   of the forward run (dE_c/dtheta: fraction of the incident Rayleigh energy radiated as c per
%   radian at theta), and the side face at 90 deg - theta.
arguments
    opts.ppw (1,1) double = 8                  % GLL intervals per lambda_S (p = 4)
    opts.corner (1,1) logical = true
    opts.vpvs (1,1) double = 4814/2962         % homogenised medium
    opts.x0 (1,1) double = -6                  % packet centre (lambda_S)
    opts.rEnd (1,1) double = 16              % Rayleigh distance from the corner at tEnd
    opts.relWidth (1,1) double = 0.25          % packet sigma_k / k_R
    opts.nBin (1,1) double = 18
    opts.Yc (1,1) double = 2.5                 % strip width / lambda_R
    opts.incident char {mustBeMember(opts.incident, {'R', 'P', 'S'})} = 'R'
    opts.theta (1,1) double = 45               % body incidence: direction it comes from (deg)
    opts.w (1,1) double = 3                    % body beam width (lambda_S)
    opts.xMin (1,1) double = NaN               % body incidence: Rayleigh strips start here (default rEnd/2)
end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', '..', 'SEMWavePropagation', 'src'));
vs = 1;  vp = opts.vpvs;  rho = 1;
nu = (vp^2 - 2)/(2*(vp^2 - 1));
mu = rho*vs^2;  lam = rho*vp^2 - 2*mu;
xi = fzero(@(x) (2 - x)^2 - 4*sqrt(1 - x/vp^2)*sqrt(1 - x), [0.5 0.99]);
cR = sqrt(xi);  qa = sqrt(1 - xi/vp^2);  sa = sqrt(1 - xi);
isBeam = opts.incident ~= 'R';
c0 = cR;  if opts.incident == 'P', c0 = vp; elseif opts.incident == 'S', c0 = vs; end
k0 = 2*pi/c0;                                                      % omega = 2 pi (lambda_S = 1)
Rp = vp*opts.rEnd/cR + 2;                                          % P front from the corner + margin
th0 = opts.theta*pi/180;  nb = [cos(th0) -sin(th0)];  tb = [sin(th0) cos(th0)];
Dist = 4 + max((1 + 3*opts.w*[cos(th0) sin(th0)])./[sin(th0) cos(th0)]);   % beam clear of the faces
if isBeam, tEnd = Dist/c0 + opts.rEnd/cR;  Rp = 1.4*opts.rEnd + 4;  else, tEnd = (-opts.x0 + opts.rEnd)/cR; end   % beam: far-face echoes miss the strips
p = 4;  h = p/opts.ppw;
if opts.corner, xb = [-ceil(max([Rp, -opts.x0 + 4, Dist + 4])) 0];  else, xb = [-(-opts.x0 + 4 + 2) opts.rEnd + 6]; end
yb = [0 ceil(max(Rp, Dist + 4))];
refBeam = isBeam && ~opts.corner;                                  % beam reference: no faces, fluence at the origin
if refBeam
    Lb = ceil(Dist + 4 + 3*opts.w);  La = ceil(3*opts.w + 8);
    xb = [-Lb La];  yb = [-La Lb];  tEnd = (Dist + 3.5*c0/(2*pi*opts.relWidth) + 2)/c0;
end
ne = round([diff(xb) diff(yb)]/h);
mesh = sem_mesh([xb; yb], ne, p);
model = sem_assemble(mesh, struct('E', 2*mu*(1 + nu), 'nu', nu, 'rho', rho), 'PSV');
op = sem_pml_operators(mesh);
fprintf('wedge SEM: vp/vs %.4f (nu %.4f), cR/vs %.4f, ppw %g (%.1f per lambda_R), %d nodes, corner %d\n', ...
        vp, nu, cR, opts.ppw, opts.ppw*cR, size(mesh.nodes, 1), opts.corner);

% exact Rayleigh packet at t = 0 and t = -dt
dt = 0.8*model.dtCritical;  nSteps = ceil(tEnd/dt);  dt = tEnd/nSteps;
X = mesh.nodes(:, 1);  Y = mesh.nodes(:, 2);
u0 = packet(0);  um = packet(-dt);
E0 = energy(u0, um, dt);

% central differences
Minv = 1./full(diag(model.M));  K = model.K;
u = u0;  uPrev = um;
if refBeam                                                         % Poynting flux along nb at the origin node
    [~, n0] = min(sum(mesh.nodes.^2, 2));  r0 = find(op.globalNodes == n0, 1);
    b1 = op.B{1}(r0, :);  b2 = op.B{2}(r0, :);  cl = find(b1 | b2);  b1 = full(b1(cl));  b2 = full(b2(cl));
    Fref = 0;
end
tic;
for n = 1:nSteps
    uNext = 2*u - uPrev - dt^2*Minv.*(K*u);
    if refBeam
        ux = u(2*cl - 1);  uy = u(2*cl);  vo = (uNext(2*n0 - [1 0]) - uPrev(2*n0 - [1 0]))/(2*dt);
        gxx = b1*ux;  gyx = b2*ux;  gxy = b1*uy;  gyy = b2*uy;     % g_ja = d u_a / d x_j
        dvo = gxx + gyy;  sxx = lam*dvo + 2*mu*gxx;  syy = lam*dvo + 2*mu*gyy;  sxy = mu*(gyx + gxy);
        Fref = Fref - dt*((sxx*vo(1) + sxy*vo(2))*nb(1) + (sxy*vo(1) + syy*vo(2))*nb(2));
    end
    uPrev = u;  u = uNext;
end
uNext = 2*u - uPrev - dt^2*Minv.*(K*u);
fprintf('  %d steps in %.0f s\n', nSteps, toc);
v = (uNext - uPrev)/(2*dt);
E1 = 0.5*v'*model.M*v + 0.5*u'*K*u;

% qp densities
U = reshape(u, 2, []).';  V = reshape(v, 2, []).';
g = cell(2, 1);  for j = 1:2, g{j} = op.B{j}*U; end              % g{j}(:, a) = d u_a / d x_j
dv = g{1}(:, 1) + g{2}(:, 2);  cu = g{1}(:, 2) - g{2}(:, 1);
exx = g{1}(:, 1);  eyy = g{2}(:, 2);  exy = (g{1}(:, 2) + g{2}(:, 1))/2;
Vq = V(op.globalNodes, :);  w = op.weights;
e = 0.5*rho*sum(Vq.^2, 2) + 0.5*(lam*dv.^2 + 2*mu*(exx.^2 + eyy.^2 + 2*exy.^2));
eP = rho*vp^2*dv.^2;  eS = rho*vs^2*cu.^2;
xq = mesh.nodes(op.globalNodes, 1);  yq = mesh.nodes(op.globalNodes, 2);
Yc = opts.Yc*cR;
inR = yq < Yc & xq < -Yc;  inT = xq > -Yc & yq > Yc;
if ~opts.corner, inR = yq < Yc;  inT = false(size(yq)); end
if isBeam
    xm = opts.xMin;  if isnan(xm), xm = opts.rEnd/2; end
    inR = yq < Yc & xq < -xm;  inT = xq > -Yc & yq > xm;
    kk = spectrum();  A2 = kk.^2.*exp(-(kk - k0).^2/opts.relWidth^2/k0^2);
    F = E0/(opts.w*sqrt(pi));                                      % fluence at the waist (start)
    gk = (1 + (Dist./(kk*opts.w^2)).^2).^-0.5;                     % 2-D Gaussian beam spreading to the corner, z_R = k w^2
    lamAvg = sum(2*pi./kk.*gk.*A2)/sum(A2);                        % <lambda> weighted by the paraxial fluence at the corner
    out = struct('incident', opts.incident, 'theta', opts.theta, 'w', opts.w, 'ppw', opts.ppw, 'lamAvg', lamAvg, ...
                 'lamF', sum(2*pi./kk.*gk.*A2)/sum(gk.*A2), 'gParax', sum(gk.*A2)/sum(A2), 'E0', E0, ...
                 'sigmaTop', sum(w(inR).*e(inR))/F, 'sigmaSide', sum(w(inT).*e(inT))/F, ...
                 'Wtop', sum(w(inR).*e(inR))/E0, 'Wside', sum(w(inT).*e(inT))/E0, 'conservation', E1/E0 - 1);
    if refBeam                                                     % measured fluence at the origin / waist fluence
        out.Fref = Fref/E0;  out.gMeas = Fref/F;
        fprintf('  %s beam w %g from %g deg, no faces: F(origin)/F(waist) %.4f (paraxial %.4f)\n', ...
                opts.incident, opts.w, opts.theta, out.gMeas, out.gParax);
        return
    end
    fprintf('  %s from %g deg: sigma top %.4f, side %.4f (lambda_S units), <lambda> %.4f, E conservation %.1e\n', ...
            opts.incident, opts.theta, out.sigmaTop, out.sigmaSide, lamAvg, out.conservation);
    return
end
body = ~inR & ~inT;
Eq = sum(w.*e);
out = struct('ppw', opts.ppw, 'corner', opts.corner, 'nu', nu, 'cR', cR, 'E0', E0, 'E1', E1, ...
             'conservation', E1/E0 - 1, 'qpVsMK', Eq/E1 - 1);
out.R = sum(w(inR).*e(inR))/E0;  out.T = sum(w(inT).*e(inT))/E0;
out.body = sum(w(body).*e(body))/E0;
out.P = sum(w(body).*eP(body))/E0;  out.S = sum(w(body).*eS(body))/E0;
out.bodySplitCheck = (out.P + out.S)/out.body - 1;
th = atan2(yq, -xq);  ed = linspace(0, pi/2, opts.nBin + 1);  out.thetaEdges = ed;
b = discretize(th(body), ed);  ok = ~isnan(b);
wbody = w(body);  ebP = eP(body);  ebS = eS(body);
out.patP = accumarray(b(ok), wbody(ok).*ebP(ok), [opts.nBin 1]).'/E0;
out.patS = accumarray(b(ok), wbody(ok).*ebS(ok), [opts.nBin 1]).'/E0;
sxx = lam*dv + 2*mu*exx;  syy = lam*dv + 2*mu*eyy;  sxy = 2*mu*exy;  % patterns by the Poynting direction
Px = -(sxx.*Vq(:, 1) + sxy.*Vq(:, 2));  Py = -(sxy.*Vq(:, 1) + syy.*Vq(:, 2));   % (no parallax from the
phi = atan2(Py(body), -Px(body));  b = discretize(phi, ed);  ok = ~isnan(b);     %  extended conversion zone)
out.patPdir = accumarray(b(ok), wbody(ok).*ebP(ok), [opts.nBin 1]).'/E0;
out.patSdir = accumarray(b(ok), wbody(ok).*ebS(ok), [opts.nBin 1]).'/E0;
out.sum = out.R + out.T + out.body;
fprintf('  E conservation %.1e, qp vs M/K %.1e\n', out.conservation, out.qpVsMK);
fprintf('  R %.4f  T %.4f  body %.4f (P %.4f, S %.4f, split check %.1e)  sum %.4f\n', ...
        out.R, out.T, out.body, out.P, out.S, out.bodySplitCheck, out.sum);

    function kk = spectrum()
        sk = opts.relWidth*k0;  kk = linspace(max(k0 - 5*sk, 0.05*k0), k0 + 5*sk, 600);
    end

    function uu = packet(t)
        sk = opts.relWidth*k0;  kk = spectrum();  dk = kk(2) - kk(1);
        ux = zeros(size(X));  uy = ux;
        if isBeam                                                  % plane P / S beam toward the corner
            p0 = -Dist*nb;  sl = (X - p0(1))*nb(1) + (Y - p0(2))*nb(2);
            m = (X - p0(1))*tb(1) + (Y - p0(2))*tb(2);
            pol = nb;  if opts.incident == 'S', pol = tb; end
            f = zeros(size(X));
            for c = 1:20:numel(kk)
                kc = kk(c:min(c + 19, end));
                f = f + real(1i*exp(1i*kc.*(sl - c0*t)))*(exp(-(kc - k0).^2/(2*sk^2))*dk).';
            end
            f = f.*exp(-m.^2/(2*opts.w^2));
            uu = reshape([pol(1)*f pol(2)*f].', [], 1);
            return
        end
        for c = 1:20:numel(kk)
            kc = kk(c:min(c + 19, end));
            A = exp(-(kc - k0).^2/(2*sk^2))*dk;  ph = exp(1i*(kc.*(X - opts.x0) - cR*kc*t));
            eq = exp(-qa*Y*kc);  es = exp(-sa*Y*kc);
            ux = ux + real((ph.*(1i*(eq - 2*qa*sa/(sa^2 + 1)*es)))*A.');
            uy = uy + real((ph.*(-qa*(eq - 2/(sa^2 + 1)*es)))*A.');
        end
        uu = reshape([ux uy].', [], 1);
    end

    function E = energy(a, aPrev, dtt)
        aNext = 2*a - aPrev - dtt^2*(1./full(diag(model.M))).*(model.K*a);
        vv = (aNext - aPrev)/(2*dtt);
        E = 0.5*vv'*model.M*vv + 0.5*a'*model.K*a;
    end
end
