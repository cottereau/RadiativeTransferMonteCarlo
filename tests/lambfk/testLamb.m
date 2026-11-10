function tests = testLamb
%TESTLAMB  Regression and physics tests for the +lamb package.
%   results = runtests('testLamb')      (from tests/lambfk)
%   Reference numbers: LambFK_Formulation.tex, Sec. 7 (independent Python evaluation).
%
%   References:
%     - Miller, G.F. & Pursey, H. (1955), On the partition of energy between elastic waves
%       in a semi-infinite solid, Proc. R. Soc. Lond. A 233, 55-69.
%       doi:10.1098/rspa.1955.0245
tests = functiontests(localfunctions);
end

function setupOnce(tc)
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));  AddRTMCLib(fullfile(fileparts(mfilename('fullpath')), '..', '..'), false);
tc.TestData.s = lamb.Medium(1, sqrt(3), 1);                  % Poisson solid
tc.TestData.M = [1 .4 -.7; .4 -.3 .5; -.7 .5 .6];            % general moment tensor
tc.TestData.F = [0.3; -0.5; 0.8];                            % oblique force
end

function testWholeSpace(tc)
% Identical half-spaces = whole space: analytic force power, S/P ratios, no poles.
s = tc.TestData.s;  I = lamb.Interface(s, s);  F = tc.TestData.F;
EP = lamb.energyPartition(I, lamb.PointSource('force', F, 0.3), 1);
Wa = norm(F)^2*(1/s.alpha^3 + 2/s.beta^3)/(24*pi*s.rho);
verifyEqual(tc, EP.W.total, Wa, 'RelTol', 1e-8);
verifyEmpty(tc, EP.W.interface);
b = EP.W.body;
verifyEqual(tc, (b.SV1 + b.SH1 + b.SV2 + b.SH2)/(b.P1 + b.P2), 2*sqrt(3)^3, 'RelTol', 1e-8);
EP = lamb.energyPartition(I, lamb.PointSource.doubleCouple(1, 30, 60, 90, 0.3), 1);
b = EP.W.body;
verifyEqual(tc, (b.SV1 + b.SH1 + b.SV2 + b.SH2)/(b.P1 + b.P2), 1.5*sqrt(3)^5, 'RelTol', 1e-8);
end

function testMillerPursey(tc)
% Vertical point force on a free Poisson half-space: 6.9 / 25.8 / 67.4 % (Miller & Pursey 1955).
EP = lamb.energyPartition(lamb.Interface(tc.TestData.s), lamb.PointSource('force', [0 0 1], 0), 1);
verifyEqual(tc, EP.pct.body.P1, 6.88, 'AbsTol', 0.01);
verifyEqual(tc, EP.pct.body.SV1 + EP.pct.body.SH1, 25.75, 'AbsTol', 0.01);
verifyEqual(tc, EP.pct.ifTotal, 67.36, 'AbsTol', 0.01);
verifyLessThan(tc, EP.check.relErr, 1e-8);
end

function testEnergyBalance(tc)
% Source work = body + interface power, for every case and source type.
s = tc.TestData.s;  M = tc.TestData.M;
EP = lamb.energyPartition(lamb.Interface(s), lamb.PointSource('moment', M, 0.5), 1);
verifyEqual(tc, [EP.pct.body.P1, EP.pct.body.SV1, EP.pct.body.SH1, EP.pct.ifTotal], ...
            [2.69 8.48 62.02 26.81], 'AbsTol', 0.01);
cases = {lamb.Interface(s, lamb.Medium(1.0, 1.5)),     lamb.PointSource('moment', M, 0.5)
         lamb.Interface(s, lamb.Medium(1.0, 0.5)),     lamb.PointSource('force', tc.TestData.F, 0.2)
         lamb.Interface(s, lamb.Medium(1.0, 1.5)),     lamb.PointSource.explosion(1, -0.5)
         lamb.Interface(s, lamb.Medium(1.0, 0.5)),     lamb.PointSource('moment', M, -0.3)
         lamb.Interface(s, lamb.Medium(2.0, 2.2, 1.2)), lamb.PointSource('moment', M, 0.4)
         lamb.Interface(s, lamb.Medium(1.1, 1.8, 1.05)),lamb.PointSource('force', tc.TestData.F, 0.1)};
for c = 1:size(cases, 1)
    EP = lamb.energyPartition(cases{c, 1}, cases{c, 2}, 1);
    verifyLessThan(tc, EP.check.relErr, 1e-7, sprintf('case %d', c));
    if ~isempty(EP.W.interface)
        verifyLessThan(tc, EP.check.ifRelErr, 1e-7, sprintf('case %d (IW)', c));
    end
end
end

function testRTUnitarity(tc)
% Plane-wave energy coefficients sum to one for propagating incidence.
s = tc.TestData.s;
ifs = {lamb.Interface(s), lamb.Interface(s, lamb.Medium(1, 1.5)), ...
       lamb.Interface(s, lamb.Medium(2, 2.2, 1.2))};
for i = 1:numel(ifs)
    inc = {'P1', 'SV1', 'SH1'};
    if strcmp(ifs{i}.kind, 'solid-fluid'), inc{end+1} = 'Pf'; end %#ok<AGROW>
    for j = 1:numel(inc)
        if strcmp(inc{j}, 'Pf'), c = ifs{i}.upper.alpha;
        elseif inc{j}(1) == 'P', c = s.alpha; else, c = s.beta; end
        p  = linspace(0.003, 0.997, 97)/c;
        rt = lamb.energy.rtCoefficients(ifs{i}, inc{j}, p, 1);
        verifyEqual(tc, rt.sum, ones(size(rt.sum)), 'AbsTol', 1e-10, ...
                    sprintf('%s, %s', ifs{i}.kind, inc{j}));
    end
end
end

function testFreeSurfaceLimit(tc)
% Weak fluid loading -> free surface: fast fluid (real Scholte -> Rayleigh) and slow fluid
% (near-real leaky Rayleigh, quasi-trapped). Both must give the Miller-Pursey partition.
s = tc.TestData.s;  src = lamb.PointSource('force', [0 0 1], 0);
for f = {lamb.Medium(1e-6, 1.5), lamb.Medium(1e-6, 1e-6)}
    EP = lamb.energyPartition(lamb.Interface(s, f{1}), src, 1);
    verifyEqual(tc, EP.pct.trapped, 67.36, 'AbsTol', 0.01);
    verifyEqual(tc, EP.pct.body.P1, 6.88, 'AbsTol', 0.01);
    verifyEqual(tc, EP.pct.body.SV1 + EP.pct.body.SH1, 25.75, 'AbsTol', 0.01);
    verifyLessThan(tc, EP.check.relErr, 1e-6);
end
end

function testDWNWholeSpace(tc)
% DWN wavefield vs the analytic whole-space Green's tensor at a complex frequency,
% for a force and (by finite differences of the Green's tensor) a moment tensor.
s = tc.TestData.s;  I = lamb.Interface(s, s);  zs = 0.3;  w = 1 + 0.2i;
rcv = struct('x', [1; -0.7; 0.3], 'y', [0.5; 2; -0.2], 'z', [1.3; 0.8; -0.9]);
X = [rcv.x rcv.y rcv.z];
F = tc.TestData.F;  M = tc.TestData.M;
U = lamb.wavefield.dwn(I, lamb.PointSource('force', F, zs), rcv, w, 'L', 250);
Ua = zeros(3, 3);  Um = zeros(3, 3);  h = 1e-4;
for r = 1:3
    Ua(r, :) = greens(s, w, X(r,:) - [0 0 zs])*F;
    for q = 1:3
        e = zeros(1, 3);  e(q) = h;
        dG = (greens(s, w, X(r,:) - [0 0 zs] - e) - greens(s, w, X(r,:) - [0 0 zs] + e))/(2*h);
        Um(r, :) = Um(r, :) + (dG*M(:, q)).';          % u_i = M_pq dG_ip/dxi_q
    end
end
verifyEqual(tc, squeeze(U), Ua, 'RelTol', 1e-6);
U = lamb.wavefield.dwn(I, lamb.PointSource('moment', M, zs), rcv, w, 'L', 250);
verifyEqual(tc, squeeze(U), Um, 'RelTol', 1e-5);
end

function testFreq2TimeRicker(tc)
% Complex-frequency synthesis of a Ricker spectrum must return the Ricker wavelet.

fp = 2;  t0 = 1;  dt = 0.005;  Nt = 1024;  T = Nt*dt;  ep = pi/T;
wc = (0:Nt/2).'*(2*pi/T) + 1i*ep;
S = lamb.PointSource.ricker(fp, t0);
[t, u] = lamb.wavefield.freq2time(S(wc), dt, ep);
a  = (pi*fp*(t - t0)).^2;
verifyEqual(tc, u, (1 - 2*a).*exp(-a), 'AbsTol', 1e-8);
end

function testLineSource2D(tc)
% Line source (lamb.energyPartition2D). Whole space: power per unit length of an in-plane line
% force F^2 w/(16 rho) (1/alpha^2 + 1/beta^2) ((w/2) Im G_zz(0), G the 2D Green's tensor,
% Im = (1/(8 rho))(1/alpha^2 + 1/beta^2)), S/P = (alpha/beta)^2. Free Poisson half-space,
% vertical line force at h = 0: P 21.94 %, S 29.13 %, Rayleigh 48.93 %.
s = tc.TestData.s;  w = 1;
EP = lamb.energyPartition2D(lamb.Interface(s, s), lamb.PointSource('force', [0 0 1], 0.3), w);
b = EP.W.body;
verifyEqual(tc, EP.W.total, w/(16*s.rho)*(1/s.alpha^2 + 1/s.beta^2), 'RelTol', 1e-8);
verifyEqual(tc, (b.SV1 + b.SV2)/(b.P1 + b.P2), 3, 'RelTol', 1e-8);
verifyEqual(tc, b.SH1 + b.SH2, 0, 'AbsTol', 1e-14);
EP = lamb.energyPartition2D(lamb.Interface(s), lamb.PointSource('force', [0 0 1], 0), w);
verifyEqual(tc, EP.pct.body.P1, 21.94, 'AbsTol', 0.01);
verifyEqual(tc, EP.pct.body.SV1, 29.13, 'AbsTol', 0.01);
verifyEqual(tc, EP.pct.ifTotal, 48.93, 'AbsTol', 0.01);
end

function G = greens(m, w, R)
% Whole-space Green's tensor G_ij = [k_b^2 d_ij g_b + d_i d_j (g_b - g_a)]/(rho w^2).
r = norm(R);  g = R(:)/r;
ka = w/m.alpha;  kb = w/m.beta;
gf  = @(k) exp(1i*k*r)/(4*pi*r);
d1  = @(k) (1i*k - 1/r)*gf(k);
d2  = @(k) ((1i*k - 1/r)^2 + 1/r^2)*gf(k);
H   = @(k) d2(k)*(g*g.') + d1(k)/r*(eye(3) - g*g.');
G   = (kb^2*eye(3)*gf(kb) + H(kb) - H(ka))/(m.rho*w^2);
end
