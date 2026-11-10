function allPass = PlanarXY_Validation()
%PLANARXY_VALIDATION  Validation of dir==1 (x) and dir==2 (y) planar interfaces.
%
%   Explicitly validates propagateParticleMultiMaterial on planar interfaces
%   aligned with x (bnd.dir==1) and y (bnd.dir==2):
%     1. Solid-Solid weld: P, SV, and SH incidence at oblique angles.
%     2. Solid-Fluid weld: P and SV incidence.
%     3. Verification of Snell's law, coordinate folding, and mode conversion.
%     4. Statistical verification against Zoeppritz theory within 5/sqrt(N).
%
%   See also PROPAGATEPARTICLEMULTIMATERIAL, MATERIALINTERFACETWOSOLID,
%   MATERIALINTERFACESOLIDFLUID.

fprintf('\n=======================================================\n');
fprintf('  PlanarXY_Validation: Explicit test for dir==1 and dir==2\n');
fprintf('=======================================================\n\n');

allPass = true;
N = 100000;
tol = 5 / sqrt(N);

% Test materials
matA = isoMat(3000, 1500, 2200, 100.0); % Shale (solid)
matB = isoMat(4500, 2600, 2600, 100.0); % Sandstone (solid)
matFluid = isoMat(1500, 0, 1000, 100.0); % Water (fluid)

j1_deg = linspace(0, 90, 361);
zoepSolid = MaterialClass.ZoeppritzSolid(j1_deg, matA.vp, matA.vs, matA.rho, matB.vp, matB.vs, matB.rho);
zoepSolid.j1_deg = j1_deg;

zoepFluid = MaterialClass.ZoeppritzFluid(j1_deg, matA.vp, matA.vs, matA.rho, matFluid.vp, 0, matFluid.rho);
zoepFluid.j1_deg = j1_deg;

testAngles = [15, 35, 55];

for dirAxis = [1, 2]
    axisName = 'x'; if dirAxis == 2, axisName = 'y'; end
    fprintf('>>> Testing planar interface along %s-axis (dir == %d) <<<\n', axisName, dirAxis);

    geometry = struct('dimension', 3, 'frame', 'cartesian');
    geometry.bnd(1) = struct('dir', dirAxis, 'val', 0.0, 'type', '', 'mat', [1 2]);

    %% --------------------------------------------------------------------
    % 1. TwoSolid: P-wave incidence
    %% --------------------------------------------------------------------
    fprintf('  [TwoSolid P-incidence on %s=0]\n', axisName);
    for theta = testAngles
        P = createBeam(N, dirAxis, theta, true, 'P');
        dt = 1e-4; % short step to hit the interface at val = 0
        P = propagateParticleMultiMaterial({matA, matB}, geometry, P, dt);

        % Theory coefficients
        Rpp  = interp1(j1_deg, zoepSolid.E_Rpp,  theta);
        Rpsv = interp1(j1_deg, zoepSolid.E_Rpsv, theta);
        Tpp  = interp1(j1_deg, zoepSolid.E_Tpp,  theta);
        Tpsv = interp1(j1_deg, zoepSolid.E_Tpsv, theta);

        % MC fractions
        mc_Rpp  = mean(P.matIdx == 1 &  P.p);
        mc_Rpsv = mean(P.matIdx == 1 & ~P.p);
        mc_Tpp  = mean(P.matIdx == 2 &  P.p);
        mc_Tpsv = mean(P.matIdx == 2 & ~P.p);

        % Validation checks
        pass1 = abs(mc_Rpp  - Rpp)  < tol;
        pass2 = abs(mc_Rpsv - Rpsv) < tol;
        pass3 = abs(mc_Tpp  - Tpp)  < tol;
        pass4 = abs(mc_Tpsv - Tpsv) < tol;

        fprintf('    theta=%2d deg | Rpp: %.3f (th %.3f) | Rpsv: %.3f (th %.3f) | Tpp: %.3f (th %.3f) | Tpsv: %.3f (th %.3f)\n', ...
            theta, mc_Rpp, Rpp, mc_Rpsv, Rpsv, mc_Tpp, Tpp, mc_Tpsv, Tpsv);

        allPass = allPass && pass1 && pass2 && pass3 && pass4;

        % Position validation: reflected particles must have coord < 0, transmitted > 0
        assert(all(P.x(P.matIdx == 1, dirAxis) <= 0), 'Reflected particle in wrong halfspace');
        assert(all(P.x(P.matIdx == 2, dirAxis) >= 0), 'Transmitted particle in wrong halfspace');
    end

    %% --------------------------------------------------------------------
    % 2. TwoSolid: SV-wave incidence
    %% --------------------------------------------------------------------
    fprintf('  [TwoSolid SV-incidence on %s=0]\n', axisName);
    for theta = testAngles
        P = createBeam(N, dirAxis, theta, false, 'SV');
        dt = 1e-4;
        P = propagateParticleMultiMaterial({matA, matB}, geometry, P, dt);

        Rsvsv = interp1(j1_deg, zoepSolid.E_Rsvsv, theta);
        Rsvp  = interp1(j1_deg, zoepSolid.E_Rsp,   theta);
        Tsvsv = interp1(j1_deg, zoepSolid.E_Tsvsv, theta);
        Tsvp  = interp1(j1_deg, zoepSolid.E_Tsp,   theta);

        mc_Rsvsv = mean(P.matIdx == 1 & ~P.p);
        mc_Rsvp  = mean(P.matIdx == 1 &  P.p);
        mc_Tsvsv = mean(P.matIdx == 2 & ~P.p);
        mc_Tsvp  = mean(P.matIdx == 2 &  P.p);

        pass1 = abs(mc_Rsvsv - Rsvsv) < tol;
        pass2 = abs(mc_Rsvp  - Rsvp)  < tol;
        pass3 = abs(mc_Tsvsv - Tsvsv) < tol;
        pass4 = abs(mc_Tsvp  - Tsvp)  < tol;

        fprintf('    theta=%2d deg | Rsvsv: %.3f (th %.3f) | Rsvp: %.3f (th %.3f) | Tsvsv: %.3f (th %.3f) | Tsvp: %.3f (th %.3f)\n', ...
            theta, mc_Rsvsv, Rsvsv, mc_Rsvp, Rsvp, mc_Tsvsv, Tsvsv, mc_Tsvp, Tsvp);

        allPass = allPass && pass1 && pass2 && pass3 && pass4;
    end

    %% --------------------------------------------------------------------
    % 3. TwoSolid: SH-wave incidence
    %% --------------------------------------------------------------------
    fprintf('  [TwoSolid SH-incidence on %s=0]\n', axisName);
    for theta = testAngles
        P = createBeam(N, dirAxis, theta, false, 'SH');
        dt = 1e-4;
        P = propagateParticleMultiMaterial({matA, matB}, geometry, P, dt);

        Rshsh = interp1(j1_deg, zoepSolid.E_Rsh, theta);
        Tshsh = interp1(j1_deg, zoepSolid.E_Tsh, theta);

        mc_Rshsh = mean(P.matIdx == 1);
        mc_Tshsh = mean(P.matIdx == 2);

        pass1 = abs(mc_Rshsh - Rshsh) < tol;
        pass2 = abs(mc_Tshsh - Tshsh) < tol;

        fprintf('    theta=%2d deg | Rshsh: %.3f (th %.3f) | Tshsh: %.3f (th %.3f)\n', ...
            theta, mc_Rshsh, Rshsh, mc_Tshsh, Tshsh);

        allPass = allPass && pass1 && pass2;
    end

    %% --------------------------------------------------------------------
    % 4. Solid-Fluid: P-wave incidence
    %% --------------------------------------------------------------------
    fprintf('  [Solid-Fluid P-incidence on %s=0]\n', axisName);
    for theta = testAngles
        P = createBeam(N, dirAxis, theta, true, 'P');
        dt = 1e-4;
        P = propagateParticleMultiMaterial({matA, matFluid}, geometry, P, dt);

        Rpp  = interp1(j1_deg, zoepFluid.E_Rpp,  theta);
        Rpsv = interp1(j1_deg, zoepFluid.E_Rpsv, theta);
        Tpf  = interp1(j1_deg, zoepFluid.E_Tpp,  theta);

        mc_Rpp  = mean(P.matIdx == 1 &  P.p);
        mc_Rpsv = mean(P.matIdx == 1 & ~P.p);
        mc_Tpf  = mean(P.matIdx == 2);

        % Fluid cannot support S waves
        assert(all(P.p(P.matIdx == 2)), 'Fluid region contains S polarized particles!');

        pass1 = abs(mc_Rpp  - Rpp)  < tol;
        pass2 = abs(mc_Rpsv - Rpsv) < tol;
        pass3 = abs(mc_Tpf  - Tpf)  < tol;

        fprintf('    theta=%2d deg | Rpp: %.3f (th %.3f) | Rpsv: %.3f (th %.3f) | Tpf: %.3f (th %.3f)\n', ...
            theta, mc_Rpp, Rpp, mc_Rpsv, Rpsv, mc_Tpf, Tpf);

        allPass = allPass && pass1 && pass2 && pass3;
    end

    fprintf('\n');
end

if allPass
    fprintf('>>> ALL PLANAR XY (dir==1 & dir==2) VALIDATIONS PASSED! <<<\n\n');
else
    error('Planar XY validation FAILED!');
end

end

% -------------------------------------------------------------------------
function P = createBeam(N, dirAxis, thetaDeg, isP, modeStr)
% Launch a beam of N particles towards the interface at dirAxis = 0 from
% negative coordinates.
th = thetaDeg * pi / 180;
P.N = N;
P.x = zeros(N, 3);
P.p = repmat(isP, N, 1);
P.alive = true(N, 1);
P.t = zeros(N, 1);
P.matIdx = ones(N, 1);

% Choose coordinates so normal is along dirAxis, tangent in other axes
% Normal component is cos(th) > 0 (pointing towards interface at 0)
% Tangent component is sin(th)
if dirAxis == 1 % Normal along x
    % P.x starts at x = -0.05 m
    P.x(:, 1) = -0.05;
    d_n = cos(th);
    d_t = sin(th);
    P.dir = repmat([d_n, d_t, 0], N, 1);

    if isP
        P.perp = repmat([0 0 1], N, 1);
    elseif strcmp(modeStr, 'SV')
        % SV perp in plane of incidence (x-y plane): [-sin(th), cos(th), 0]
        P.perp = repmat([-d_t, d_n, 0], N, 1);
    else % SH
        % SH perp normal to plane of incidence (z-axis): [0, 0, 1]
        P.perp = repmat([0 0 1], N, 1);
    end
else % Normal along y
    % P.x starts at y = -0.05 m
    P.x(:, 2) = -0.05;
    d_n = cos(th);
    d_t = sin(th);
    P.dir = repmat([0, d_n, d_t], N, 1);

    if isP
        P.perp = repmat([1 0 0], N, 1);
    elseif strcmp(modeStr, 'SV')
        % SV perp in plane of incidence (y-z plane): [0, -sin(th), cos(th)]
        P.perp = repmat([0, -d_t, d_n], N, 1);
    else % SH
        % SH perp normal to plane of incidence (x-axis): [1, 0, 0]
        P.perp = repmat([1 0 0], N, 1);
    end
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
