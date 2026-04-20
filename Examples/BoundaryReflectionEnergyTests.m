%
% Validation of elastic reflecting-boundary treatment within the code.
%
% The test uses a single planar boundary z = 0 and incident plane-wave
% packets of type P, pure SV, and pure SH at several incidence angles.
%
% No volume scattering is considered. Only the boundary reflection and
% mode conversions are tested.
%
% Output:
%   1. Time histories of E_P, E_SV, E_SH, E_total
%   2. A printed table comparing measured final energies with:
%      - code expectation
%      - reflection-only normalized expectation
%
% Important:
%   The usual obs.energy separates only P and S. Here we additionally
%   separate SV and SH by projecting the S polarization onto the
%   SV/SH basis relative to the reflecting plane.

close all;
clc;
rng(1);

% -------------------------------------------------------------------------
% Numerical parameters
% -------------------------------------------------------------------------
Nparticles = 2e5;

vp  = 2500;
vs  = vp / sqrt(3);
rho = 2000;

zCenter  = 80;      % center height of the incident wave packet
aperture = 20;      % finite aperture of the plane packet
nTime    = 260;     % number of output times

% Incidence angles, in degrees, measured from the normal to the boundary.
% For SV -> P conversion, stay below the critical angle of arcsin(vs/vp) ~
% 35.26 degress for the tests
pAngles  = [15 30 45 60];
svAngles = [10 20 30];
shAngles = [15 45 60];

% -------------------------------------------------------------------------
% Geometry: one reflecting boundary z = 0
% -------------------------------------------------------------------------
geometry = struct('dimension', 3, 'frame', 'cartesian');
geometry.bnd(1) = struct('dir', 3, 'val', 0, 'type', 'reflective');

% -------------------------------------------------------------------------
% Material: 3D homogeneous elastic medium, no volume scattering
% -------------------------------------------------------------------------
% Mean free times of P & S waves are both Inf => no scattering
mat = makeNoScatteringElasticMaterial(vp, vs, rho);

% Theoretical coefficients from the code's Zoeppritz routine
[outZ, ~] = MaterialClass.Zoeppritz(mat);

thetaCritSVtoP = asind(vs / vp);

fprintf('\n============================================================\n');
fprintf('Elastic reflecting-boundary energy test\n');
fprintf('Single reflecting boundary: z = 0\n');
fprintf('vp = %.6g, vs = %.6g, rho = %.6g\n', vp, vs, rho);
fprintf('SV -> P critical angle = %.4f deg\n', thetaCritSVtoP);
fprintf('Number of particles = %d\n', Nparticles);
fprintf('============================================================\n\n');

% -------------------------------------------------------------------------
% Build list of cases
% -------------------------------------------------------------------------
caseMode  = {};
caseAngle = [];

for a = pAngles
    caseMode{end+1} = 'P'; %#ok<AGROW>
    caseAngle(end+1) = a; %#ok<AGROW>
end

for a = svAngles
    caseMode{end+1} = 'SV'; %#ok<AGROW>
    caseAngle(end+1) = a; %#ok<AGROW>
end

for a = shAngles
    caseMode{end+1} = 'SH'; %#ok<AGROW>
    caseAngle(end+1) = a; %#ok<AGROW>
end

nCases = numel(caseAngle);

% Storage for summary table
rows = cell(nCases, 17);

% -------------------------------------------------------------------------
% Plot setup
% -------------------------------------------------------------------------
figure('Name', 'Boundary reflection energy histories', ...
       'Color', 'w', ...
       'Position', [100 60 1300 900]);

nCols = 2;
nRows = ceil(nCases / nCols);
tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

% -------------------------------------------------------------------------
% Main loop over cases
% -------------------------------------------------------------------------
for iCase = 1:nCases

    modeIn   = caseMode{iCase};
    thetaDeg = caseAngle(iCase);

    fprintf('Running case %2d/%2d: incident %s, theta = %.2f deg\n', ...
            iCase, nCases, modeIn, thetaDeg);

    % Create a pure incident plane-wave packet
    [P, ~, tHitCenter, tHitMax] = createIncidentPlanePacket( ...
        Nparticles, modeIn, thetaDeg, zCenter, aperture, mat);

    % Output time vector. We run until all particles in the finite aperture
    % have crossed and reflected.
    tEnd = 1.25 * tHitMax;
    tVec = linspace(0, tEnd, nTime);

    energies = zeros(nTime, 4);
    % columns:
    % 1: E_P
    % 2: E_SV
    % 3: E_SH
    % 4: E_total

    nPlane = [0 0 1];

    % Energy at t = 0
    energies(1,:) = computeModeEnergies(P, nPlane);

    % Propagation and energy recording
    for it = 2:nTime
        P = propagateParticleSmallDt(mat, geometry, P, tVec(it));
        energies(it,:) = computeModeEnergies(P, nPlane);
    end

    % Measured final energy after all particles have reflected
    idxFinal = tVec > 1.10 * tHitMax;
    measFinal = mean(energies(idxFinal,:), 1);

    % Theoretical expectations
    [expCode, expReflOnly, rawCoeff] = expectedBoundaryEnergies( ...
        modeIn, thetaDeg, mat, outZ);

    % Error versus the actual implementation
    errCode = max(abs(measFinal(1:3) - expCode(1:3)));

    % ---------------------------------------------------------------------
    % Plot energy history
    % ---------------------------------------------------------------------
    nexttile;

    plot(tVec / tHitCenter, energies(:,1), 'r-', 'LineWidth', 1.4);
    hold on;
    plot(tVec / tHitCenter, energies(:,2), 'b-', 'LineWidth', 1.4);
    plot(tVec / tHitCenter, energies(:,3), 'g-', 'LineWidth', 1.4);
    plot(tVec / tHitCenter, energies(:,4), 'k--', 'LineWidth', 1.2);

    xline(1, 'k:', 'LineWidth', 1.0);

    xlabel('t / t_{hit,center}');
    ylabel('energy fraction');
    title(sprintf('%s incidence, \\theta = %.0f^\\circ', modeIn, thetaDeg));

    ylim([-0.05 1.05]);
    grid on;
    box on;

    if iCase == 1
        legend('P', 'SV', 'SH', 'total', ...
               'Location', 'best', ...
               'FontSize', 8);
    end

    % ---------------------------------------------------------------------
    % Store table row
    % ---------------------------------------------------------------------
    rows(iCase,:) = { ...
        modeIn, thetaDeg, ...
        measFinal(1), measFinal(2), measFinal(3), measFinal(4), ...
        expCode(1), expCode(2), expCode(3), expCode(4), ...
        expReflOnly(1), expReflOnly(2), expReflOnly(3), expReflOnly(4), ...
        errCode, rawCoeff(1), rawCoeff(2)};

end

% -------------------------------------------------------------------------
% Summary table
% -------------------------------------------------------------------------
summaryTable = cell2table(rows, ...
    'VariableNames', { ...
    'IncidentMode', 'ThetaDeg', ...
    'Measured_EP', 'Measured_ESV', 'Measured_ESH', 'Measured_Etot', ...
    'Expected_EP', 'Expected_ESV', 'Expected_ESH', 'Expected_Etot', ...
    'ReflOnlyExp_EP', 'ReflOnlyExp_ESV', 'ReflOnlyExp_ESH', 'ReflOnlyExp_Etot', ...
    'MaxAbsErr_vs_Expected', ...
    'RawCoeff_1', 'RawCoeff_2'});

fprintf('\n============================================================\n');
fprintf('Summary table\n');
fprintf('RawCoeff_1 and RawCoeff_2 mean:\n');
fprintf('  P  incidence: RawCoeff_1 = E_Rpp,   RawCoeff_2 = E_Rpsv\n');
fprintf('  SV incidence: RawCoeff_1 = E_Rsvsv, RawCoeff_2 = E_Rsp\n');
fprintf('  SH incidence: RawCoeff_1 = E_Rsh,   RawCoeff_2 = 0\n');
fprintf('============================================================\n\n');

disp(summaryTable);

% =========================================================================
% Create a 3D homogeneous elastic material without volume scattering
% =========================================================================
function mat = makeNoScatteringElasticMaterial(vp, vs, rho)

mat = MaterialClass();

mat.d = 3;
mat.acoustics = false;

mat.vp = vp;
mat.vs = vs;
mat.rho = rho;

mat.timeSteps = 0;

% No volume scattering.
mat.meanFreeTime = [Inf; Inf];

% These fields are needed by propagateParticleSmallDt, but since
% meanFreeTime = Inf, no volume scattering actually occurs.
mat.P2P = 1;
mat.S2S = 1;

zeroInv = @(u) zeros(size(u));

mat.invcdf = cell(2,2);
mat.invcdf{1,1} = zeroInv;
mat.invcdf{1,2} = zeroInv;
mat.invcdf{2,1} = zeroInv;
mat.invcdf{2,2} = zeroInv;

end

% =========================================================================
% Create a pure incident plane-wave packet
% =========================================================================
function [P, vInc, tHitCenter, tHitMax] = createIncidentPlanePacket( ...
    N, modeIn, thetaDeg, zCenter, aperture, mat)

% Boundary normal
nPlane = [0 0 1];

% Incident direction. The wave travels downward toward z = 0.
dIn = [sind(thetaDeg), 0, -cosd(thetaDeg)];
dIn = dIn / norm(dIn);

% Build a finite plane aperture perpendicular to dIn.
uv = null(dIn);
u = uv(:,1).';
v = uv(:,2).';

c1 = (rand(N,1) - 0.5) * aperture;
c2 = (rand(N,1) - 0.5) * aperture;

xCenter = [0 0 zCenter];
x = xCenter + c1 .* u + c2 .* v;

if any(x(:,3) <= 0)
    error(['Some particles are initialized below the boundary. ', ...
           'Increase zCenter or reduce aperture.']);
end

dir = repmat(dIn, N, 1);

% SV/SH basis for the incident direction
[eSV, eSH] = svshBasis(dir, nPlane);

switch upper(modeIn)

    case 'P'
        p = true(N,1);
        perp = eSH;     % irrelevant for P, but must be orthogonal
        vInc = mat.vp;

    case 'SV'
        p = false(N,1);
        perp = eSV;     % pure SV
        vInc = mat.vs;

    case 'SH'
        p = false(N,1);
        perp = eSH;     % pure SH
        vInc = mat.vs;

    otherwise
        error('Unknown incident mode. Use P, SV, or SH.');
end

P = struct('d', 3, ...
           'N', N, ...
           'x', x, ...
           'dir', dir, ...
           'perp', perp, ...
           'p', p, ...
           't', zeros(N,1));

% Hit times for each particle.
% Since dir_z < 0:
% z(t) = z0 + vInc * dir_z * t.
hitTimes = x(:,3) ./ (-vInc * dir(:,3));

tHitCenter = zCenter / (vInc * cosd(thetaDeg));
tHitMax    = max(hitTimes);

end

% =========================================================================
% Compute P, SV, SH energy fractions from particles
% =========================================================================
function E = computeModeEnergies(P, nPlane)

EP = mean(P.p);

isS = ~P.p;

ESV = 0;
ESH = 0;

if any(isS)
    dS = P.dir(isS,:);
    pS = P.perp(isS,:);

    [eSV, eSH] = svshBasis(dS, nPlane);

    aSV = dot(pS, eSV, 2);
    aSH = dot(pS, eSH, 2);

    % Sum over all particles, not only S particles, so that ESV and ESH
    % are global energy fractions.
    ESV = sum(aSV.^2) / P.N;
    ESH = sum(aSH.^2) / P.N;
end

Etot = EP + ESV + ESH;

E = [EP, ESV, ESH, Etot];

end

% =========================================================================
% SV/SH basis relative to the boundary normal
% =========================================================================
function [eSV, eSH] = svshBasis(d, nPlane)

nMat = repmat(nPlane(:).', size(d,1), 1);

% SH direction: perpendicular to the incidence plane spanned by d and n.
eSH = cross(d, nMat, 2);
normSH = vecnorm(eSH, 2, 2);

if any(normSH < 1e-12)
    error(['SV/SH basis is singular near normal incidence. ', ...
           'Use a nonzero incidence angle for SV/SH tests.']);
end

eSH = eSH ./ normSH;

% SV direction: in the incidence plane and orthogonal to d.
eSV = cross(eSH, d, 2);
eSV = eSV ./ vecnorm(eSV, 2, 2);

end

% =========================================================================
% Expected final energy fractions
% =========================================================================
function [expCode, expReflOnly, rawCoeff] = expectedBoundaryEnergies( ...
    modeIn, thetaDeg, mat, outZ)

switch upper(modeIn)

    case 'P'
        E_Rpp  = interp1(outZ.j1_deg, outZ.E_Rpp,  thetaDeg, 'linear');
        E_Rpsv = interp1(outZ.j1_deg, outZ.E_Rpsv, thetaDeg, 'linear');

        % What the code effectively uses:
        % P -> P with probability E_Rpp
        % P -> SV with probability 1 - E_Rpp
        expCode = [E_Rpp, 1 - E_Rpp, 0, 1];

        % Reflection-only normalized version:
        Rsum = E_Rpp + E_Rpsv;
        expReflOnly = [E_Rpp/Rsum, E_Rpsv/Rsum, 0, 1];

        rawCoeff = [E_Rpp, E_Rpsv];

    case 'SV'
        E_Rsvsv = interp1(outZ.j1_deg, outZ.E_Rsvsv, thetaDeg, 'linear');
        E_Rsp   = interp1(outZ.j1_deg, outZ.E_Rsp,   thetaDeg, 'linear');

        thetaCrit = asind(mat.vs / mat.vp);

        if thetaDeg < thetaCrit
            % What the code effectively uses:
            % SV -> SV with probability E_Rsvsv
            % SV -> P with probability 1 - E_Rsvsv
            expCode = [1 - E_Rsvsv, E_Rsvsv, 0, 1];

            % Reflection-only normalized version:
            Rsum = E_Rsvsv + E_Rsp;
            expReflOnly = [E_Rsp/Rsum, E_Rsvsv/Rsum, 0, 1];
        else
            % Above the SV -> P critical angle, the code attempts SV -> P
            % but then rejects it as invalid and falls back to SV reflection.
            expCode = [0, 1, 0, 1];
            expReflOnly = [0, 1, 0, 1];
        end

        rawCoeff = [E_Rsvsv, E_Rsp];

    case 'SH'
        E_Rsh = interp1(outZ.j1_deg, outZ.E_Rsh, thetaDeg, 'linear');

        expCode = [0, 0, 1, 1];
        expReflOnly = [0, 0, 1, 1];

        rawCoeff = [E_Rsh, 0];

    otherwise
        error('Unknown incident mode. Use P, SV, or SH.');
end

end