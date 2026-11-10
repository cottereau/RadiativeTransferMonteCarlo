function allPass = ValidationPlaneWaveTwoLayer()
%VALIDATIONPLANEWAVETWOLAYER  End-to-end plane wave validation through RTE pipeline.
%
%   Launches an oblique plane wave (source.type = 'plane') at an interface
%   between two elastic media, and checks from the spatial energy distribution
%   that:
%     1. The reflected and transmitted wavefronts propagate at the exact
%        Snell's law angles (specular reflection and Snell refraction).
%     2. The observed integrated energy fractions in the reflected and
%        transmitted packets match the Zoeppritz coefficients.
%
%   See also PROPAGATEPARTICLEMULTIMATERIAL, ZOEPPRITZSOLID, INITIALIZEPARTICLE.

fprintf('\n=================================================================\n');
fprintf('  ValidationPlaneWaveTwoLayer: End-to-End Plane Wave RTE Benchmark\n');
fprintf('=================================================================\n\n');

allPass = true;

% Material 1 (z < 0): Shale
mat1 = isoMat(3000, 1500, 2200, 100.0);
% Material 2 (z > 0): Sandstone
mat2 = isoMat(4500, 2600, 2600, 100.0);
matArray = {mat1, mat2};

% Interface at z = 0
geometry.dimension = 3;
geometry.frame = 'cartesian';
geometry.bnd(1) = struct('dir', 3, 'val', 0.0, 'type', '', 'mat', [1 2]);

% Incidence angle: theta = 25 degrees
theta_deg = 25;
theta = theta_deg * pi / 180;
inc_dir = [sin(theta), 0, cos(theta)];

% Source: plane wave launched from z0 = -0.06 m towards interface at z = 0
N = 100000;
source.type = 'plane';
source.direction = inc_dir;
source.position = [0, 0, -0.06];
source.extent = [0.03, 0.03];
source.lambda = 1e-3;
source.numberParticles = N;
source.polarization = 'P';

% Theoretical Zoeppritz coefficients at 25 degrees
j1_deg = linspace(0, 90, 361);
zoep = MaterialClass.ZoeppritzSolid(j1_deg, mat1.vp, mat1.vs, mat1.rho, mat2.vp, mat2.vs, mat2.rho);
Rpp_th  = interp1(j1_deg, zoep.E_Rpp,  theta_deg);
Rpsv_th = interp1(j1_deg, zoep.E_Rpsv, theta_deg);
Tpp_th  = interp1(j1_deg, zoep.E_Tpp,  theta_deg);
Tpsv_th = interp1(j1_deg, zoep.E_Tpsv, theta_deg);

% Snell angles
snell_Rpp_deg  = theta_deg;
snell_Rpsv_deg = asind(mat1.vs / mat1.vp * sind(theta_deg));
snell_Tpp_deg  = asind(mat2.vp / mat1.vp * sind(theta_deg));
snell_Tpsv_deg = asind(mat2.vs / mat1.vp * sind(theta_deg));

fprintf('Theoretical Zoeppritz Coefficients at theta = %d deg:\n', theta_deg);
fprintf('  Rpp:  %.4f (angle: %.1f deg)\n', Rpp_th,  snell_Rpp_deg);
fprintf('  Rpsv: %.4f (angle: %.1f deg)\n', Rpsv_th, snell_Rpsv_deg);
fprintf('  Tpp:  %.4f (angle: %.1f deg)\n', Tpp_th,  snell_Tpp_deg);
fprintf('  Tpsv: %.4f (angle: %.1f deg)\n\n', Tpsv_th, snell_Tpsv_deg);

%% 1. Initialize Particles
P = initializeParticle(N, 3, false, source);
P.matIdx = ones(N, 1);

%% 2. Propagate through the interface across dt sub-steps
% Total time: travel to interface and emerge on both sides
T_total = 4.0e-5; % s (P-wave travels ~0.12 m, comfortably crossing z=0)
dt = 2.0e-6;      % 20 sub-steps
Nt = round(T_total / dt);

fprintf('Propagating %d particles over %d time steps (dt = %.1e s)...\n', N, Nt, dt);
for it = 1:Nt
    P = propagateParticleMultiMaterial(matArray, geometry, P, it * dt);
end

%% 3. Measure Observed Packet Fractions
% Particles with z < 0 are reflected, z > 0 are transmitted
refl_P  = P.matIdx == 1 &  P.p;
refl_SV = P.matIdx == 1 & ~P.p;
tran_P  = P.matIdx == 2 &  P.p;
tran_SV = P.matIdx == 2 & ~P.p;

Rpp_mc  = mean(refl_P);
Rpsv_mc = mean(refl_SV);
Tpp_mc  = mean(tran_P);
Tpsv_mc = mean(tran_SV);

fprintf('\nObserved Monte Carlo Energy Fractions:\n');
fprintf('  Rpp:  MC = %.4f | Theory = %.4f | diff = %.4f\n', Rpp_mc,  Rpp_th,  abs(Rpp_mc  - Rpp_th));
fprintf('  Rpsv: MC = %.4f | Theory = %.4f | diff = %.4f\n', Rpsv_mc, Rpsv_th, abs(Rpsv_mc - Rpsv_th));
fprintf('  Tpp:  MC = %.4f | Theory = %.4f | diff = %.4f\n', Tpp_mc,  Tpp_th,  abs(Tpp_mc  - Tpp_th));
fprintf('  Tpsv: MC = %.4f | Theory = %.4f | diff = %.4f\n', Tpsv_mc, Tpsv_th, abs(Tpsv_mc - Tpsv_th));

tol = 5 / sqrt(N);
pass_amp = (abs(Rpp_mc  - Rpp_th)  < tol) && ...
           (abs(Rpsv_mc - Rpsv_th) < tol) && ...
           (abs(Tpp_mc  - Tpp_th)  < tol) && ...
           (abs(Tpsv_mc - Tpsv_th) < tol);

assert(pass_amp, 'Observed energy fractions differ from Zoeppritz theory beyond statistical tolerance');
fprintf('  -> Energy fractions PASS (within 5/sqrt(N) = %.4f)\n\n', tol);

%% 4. Measure Observed Propagation Angles (Snell's Law)
% Angle in x-z plane: atan2d(dx, |dz|)
calcAngle = @(idx) median(atan2d(P.dir(idx, 1), abs(P.dir(idx, 3))));

mc_ang_Rpp  = calcAngle(refl_P);
mc_ang_Rpsv = calcAngle(refl_SV);
mc_ang_Tpp  = calcAngle(tran_P);
mc_ang_Tpsv = calcAngle(tran_SV);

fprintf('Observed Propagation Angles vs Snell Theory:\n');
fprintf('  Reflected P:   Obs = %.2f deg | Snell = %.2f deg | diff = %.2f deg\n', ...
    mc_ang_Rpp,  snell_Rpp_deg,  abs(mc_ang_Rpp  - snell_Rpp_deg));
fprintf('  Reflected SV:  Obs = %.2f deg | Snell = %.2f deg | diff = %.2f deg\n', ...
    mc_ang_Rpsv, snell_Rpsv_deg, abs(mc_ang_Rpsv - snell_Rpsv_deg));
fprintf('  Transmitted P: Obs = %.2f deg | Snell = %.2f deg | diff = %.2f deg\n', ...
    mc_ang_Tpp,  snell_Tpp_deg,  abs(mc_ang_Tpp  - snell_Tpp_deg));
fprintf('  Transmitted SV:Obs = %.2f deg | Snell = %.2f deg | diff = %.2f deg\n', ...
    mc_ang_Tpsv, snell_Tpsv_deg, abs(mc_ang_Tpsv - snell_Tpsv_deg));

ang_tol = 0.5; % deg
pass_ang = (abs(mc_ang_Rpp  - snell_Rpp_deg)  < ang_tol) && ...
           (abs(mc_ang_Rpsv - snell_Rpsv_deg) < ang_tol) && ...
           (abs(mc_ang_Tpp  - snell_Tpp_deg)  < ang_tol) && ...
           (abs(mc_ang_Tpsv - snell_Tpsv_deg) < ang_tol);

assert(pass_ang, 'Observed propagation angles deviate from Snell''s law');
fprintf('  -> Snell propagation angles PASS (within %.1f deg)\n\n', ang_tol);

allPass = pass_amp && pass_ang;
if allPass
    fprintf('>>> END-TO-END PLANE WAVE RTE VALIDATION PASSED! <<<\n\n');
end

end

% -------------------------------------------------------------------------
function mat = isoMat(vp, vs, rho, mfp)
mat = MaterialClass;
mat.acoustics = false;
mat.vp = vp;
mat.vs = vs;
mat.rho = rho;
c = 1 / (4 * pi * mfp);
mat.sigma = {@(th)c*ones(size(th)), @(th)c*ones(size(th)); ...
             @(th)c*ones(size(th)), @(th)c*ones(size(th))};
mat = MaterialClass.prepareSigma(mat, 3);
end
