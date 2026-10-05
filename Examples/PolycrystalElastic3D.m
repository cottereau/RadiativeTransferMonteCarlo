close all
clearvars
clc

%POLYCRYSTALELASTIC3D Elastic RTE simulation in a cylindrical polycrystal.
%
% The example represents a through-transmission-style experiment in a
% cylindrical alpha-iron sample. A compact point source near the lower
% face emits P waves into the upper hemisphere. Scattering inside the grains
% produces P-P, P-S, S-P, and S-S events. The sample surfaces are treated
% as traction-free reflective boundaries.
%
% The script produces:
%   1. P, S, and total energy-density histories at one receiver;
%   2. P, S, and total energy-density maps at one selected time;
%   3. GIF and AVI movies of the three energy-density fields versus time.
%
% All material and geometrical quantities below use SI units.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
rng(1); % Reproducible Monte Carlo realization

%% 1. Numerical and display controls
% Use 5e5 particles first. Increasing this to 5e6 reduces Monte Carlo
% noise but is not required to remove pixelation from the displayed movie.
makeMovie = true;
numberParticles = 2e6;
radialBinSize = 0.025e-3;
axialBinSize = 0.025e-3;
outputTimeStep = 0.005e-6;
endTime = 0.70e-6;
interpolationFactor = 4;

%% 2. Polycrystal inputs
% Built-in names and symbols are case-insensitive. For example, replace
% 'alphaIron' by 'Al', 'aluminum', 'ALUMINIUM', 'Cu', or 'alphaTitanium'.
singleCrystalInput = 'Li';

% To provide properties manually instead, comment the line above and
% uncomment this block. These values were used by the earlier example.
% c11 = 226e9; c12 = 140e9; c44 = 116e9; rho = 7870;
% Cloc = [c11 c12 c12 0 0 0; c12 c11 c12 0 0 0; ...
%         c12 c12 c11 0 0 0; 0 0 0 c44 0 0; ...
%         0 0 0 0 c44 0; 0 0 0 0 0 c44];
% singleCrystalInput = struct( ...
%     'name','Custom alpha iron','symmetry','cubic', ...
%     'rho',rho,'C',Cloc);

% At 47 MHz, the P and S wavelengths bracket the average grain diameter defined below.
frequency = 47e6; % Hz

% average grain diameter
averageGrainDiameter = 100e-6;

% Exponential TPCF: eta(r) = exp(-r/a). Choosing a = D/(2*pi) places the
% spectral transition a*q = 1 near q = 2*pi/D. This is an explicit model
% choice, not a universal conversion between diameter and TPCF length.
correlationLength = averageGrainDiameter/(2*pi);
tpcfModel = 'exp';
tpcfParameter = correlationLength;

% To include a grain-size distribution, comment the two lines above and
% uncomment the following block. Available distribution types are
% 'lognormal', 'gamma', 'weibull', and 'truncatedNormal'. The same structure
% can be used with the 'Arguelles', 'Sha', and 'ShengKhazaie' TPCF models.
% grainSizeDistribution = struct( ...
%     'type','lognormal', ...
%     'meanDiameter',averageGrainDiameter, ...
%     'standardDeviation',0.35*averageGrainDiameter);
% tpcfModel = 'ShengKhazaie';
% tpcfParameter = grainSizeDistribution;

% Optional intrinsic attenuation. Use [QP QS] = [Inf Inf] to isolate
% scattering attenuation, as done here.
qualityFactors = [Inf Inf];

%% 3. Gometry : A cylindrical sample
sampleRadius = 1.5e-3;
sampleHeight = 3.0e-3;

geometry = struct('dimension', 3, 'frame', 'cylindrical');

geometry.bnd(1) = struct('dir', 3, 'val', 0, 'type', 'reflective');
geometry.bnd(2) = struct('dir', 3, 'val', sampleHeight, 'type', 'reflective');
geometry.bnd(3) = struct('dir', 4, 'val', sampleRadius, 'type', 'reflective');

%% 4. Compact P-wave point source
% The source is centered on the lower face. The 'upper' option keeps both
% its finite spatial packet and its directions inside the upper half-space.
% source.lambda is the spatial width of that initial packet; it is not the
% physical P wavelength, which is material.vp/frequency.
source = struct('numberParticles', numberParticles, ...
                'type', 'point', ...
                'position', [0 0 0], ...
                'direction', 'upper', ...
                'polarization', 'P', ...
                'lambda', 0.035e-3);

%% 5. Construct the polycrystal material
% Rows are incident modes and columns are scattered modes:
%
%   material.sigma = {P->P, P->S; S->P, S->S}.
%
material = MaterialClass.polycrystal3D( ...
    geometry,frequency,singleCrystalInput,tpcfModel,tpcfParameter);
material.Q = qualityFactors;
material.timeSteps = 0; % Required when reflective boundaries are present

crystal = material.singleCrystal;
velocities = [material.vp material.vs];

fprintf('\nCylindrical 3D polycrystal case\n');
fprintf('  material              = %s (%s)\n', crystal.name, crystal.symmetry);
fprintf('  sample diameter       = %.1f mm\n', 2*sampleRadius*1e3);
fprintf('  sample length         = %.1f mm\n', sampleHeight*1e3);
fprintf('  frequency             = %.2f MHz\n', frequency*1e-6);
fprintf('  TPCF model            = %s\n', material.TPCF.model);
if strcmp(material.TPCF.model,'exponential')
    fprintf('  exponential length a  = %.1f micrometres\n', ...
        material.TPCF.parameters.correlationLength*1e6);
else
    fprintf('  grain-size std. dev.  = %.1f micrometres\n', ...
        material.TPCF.parameters.standardDeviation*1e6);
end
fprintf('  average grain diameter= %.1f micrometres\n', averageGrainDiameter*1e6);
fprintf('  Vp, Vs                = %.1f, %.1f m/s\n', velocities);
fprintf('  P, S wavelength       = %.1f, %.1f micrometres\n', ...
    velocities/frequency*1e6);
fprintf('  Monte Carlo particles = %d\n\n', numberParticles);

%% 6. Observation grid in the cylindrical r-z plane
% The direct axial P front reaches the top at about 0.48 microseconds. The
% final frames show only the start of the first top reflection.
observation = struct( ...
    'x', 0:radialBinSize:sampleRadius, ...
    'y', [-pi pi], ...
    'z', 0:axialBinSize:sampleHeight, ...
    'directions', [0 pi], ...
    'time', 0:outputTimeStep:endTime);

%% 7. Run the radiative-transfer Monte Carlo solver
obs = radiativeTransfer(geometry, source, material, observation);

fprintf('  P, S total cross-section = %.4g, %.4g 1/s\n', ...
    sum(material.Sigma,2));
fprintf('  P, S mean free path   = %.2f, %.2f mm\n', ...
    material.meanFreePath*1e3);

PenergyDensity = squeeze(obs.energyDensity(:,:,:,1));
SenergyDensity = squeeze(obs.energyDensity(:,:,:,2));
totalEnergyDensity = PenergyDensity + SenergyDensity;

maximumEnergyDensity = max(totalEnergyDensity(:));
energyFloor = max(maximumEnergyDensity*1e-6, realmin);
logColorLimits = log10(maximumEnergyDensity) + [-5 0];

%% 8. Plot P, S, and total energy-density histories at one receiver
receiverRadius = 0;
receiverDepth = 0.80*sampleHeight;
[~, radialIndex] = min(abs(obs.x - receiverRadius));
[~, axialIndex] = min(abs(obs.z - receiverDepth));

receiverP = squeeze(PenergyDensity(radialIndex,axialIndex,:));
receiverS = squeeze(SenergyDensity(radialIndex,axialIndex,:));
receiverTotal = receiverP + receiverS;

figure('Color', 'w', 'Name', 'Polycrystal receiver histories');
semilogy(obs.t*1e6, max(receiverP,energyFloor), 'LineWidth', 1.5);
hold on
semilogy(obs.t*1e6, max(receiverS,energyFloor), 'LineWidth', 1.5);
semilogy(obs.t*1e6, max(receiverTotal,energyFloor), 'k-', 'LineWidth', 1.5);
grid on
box on
xlabel('Time [microseconds]');
ylabel('Normalized energy density [m^{-3}]');
title(sprintf('Receiver at r = %.2f mm, z = %.2f mm', ...
    obs.x(radialIndex)*1e3, obs.z(axialIndex)*1e3));
legend('P', 'S', 'Total', 'Location', 'best');

%% 9. Plot P, S, and total energy-density maps at one time
snapshotTime = 0.45e-6;
[~, snapshotIndex] = min(abs(obs.t - snapshotTime));

% Mirror the axisymmetric r-z result to display the full cylinder width.
diameter = [-obs.x(end:-1:1) obs.x]*1e3;
Pmap = [PenergyDensity(end:-1:1,:,snapshotIndex); ...
        PenergyDensity(:,:,snapshotIndex)];
Smap = [SenergyDensity(end:-1:1,:,snapshotIndex); ...
        SenergyDensity(:,:,snapshotIndex)];
totalMap = Pmap + Smap;

figure('Color', 'w', 'Name', 'Polycrystal energy-density snapshot');
layout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
energyMaps = {Pmap, Smap, totalMap};
panelTitles = {'P energy density', 'S energy density', 'Total energy density'};

% Refine only the displayed maps. This does not modify obs.energyDensity.
displayDiameter = linspace(diameter(1),diameter(end), ...
    interpolationFactor*(numel(diameter)-1)+1);
displayDepth = linspace(obs.z(1)*1e3,obs.z(end)*1e3, ...
    interpolationFactor*(numel(obs.z)-1)+1);
[diameterQuery,depthQuery] = ndgrid(displayDiameter,displayDepth);

for panelIndex = 1:3
    interpolant = griddedInterpolant({diameter,obs.z*1e3}, ...
        energyMaps{panelIndex},'linear','nearest');
    energyMaps{panelIndex} = interpolant(diameterQuery,depthQuery);
end

for panelIndex = 1:3
    ax = nexttile(layout);
    imagesc(ax, displayDiameter, displayDepth, ...
        log10(max(energyMaps{panelIndex},energyFloor))');
    set(ax,'YDir','normal');
    axis(ax,'equal');
    xlim(ax,[-sampleRadius sampleRadius]*1e3);
    ylim(ax,[0 sampleHeight]*1e3);
    clim(ax,logColorLimits);
    colormap(ax,'turbo');
    colorbarHandle = colorbar(ax);
    colorbarHandle.Label.String = ...
        'log_{10}(normalized energy density [m^{-3}])';
    box(ax,'on');
    ylabel(ax,'z [mm]');
    xlabel(ax,'diameter coordinate [mm]');
    title(ax,sprintf('%s, t = %.2f microseconds', ...
        panelTitles{panelIndex},obs.t(snapshotIndex)*1e6));
end

%% 10. Write GIF and AVI movies
if makeMovie
    movieFile = fullfile(repoRoot,'PolycrystalElastic3D_log');

    plotting = struct( ...
        'movieTotalEnergy', true, ...
        'energyScale', 'log10', ...
        'logFloor', energyFloor, ...
        'clim', [7 9], ...
        'filename', movieFile, ...
        'movieFormat', ["gif","avi"], ...
        'frameRate', 12.5, ...
        'figureSize', [1920 1080], ...
        'timeScale', 1e6, ...
        'timeUnit', 'microseconds', ...
        'coordinateScale', 1e3, ...
        'coordinateUnit', 'mm', ...
        'visible', 'on', ...
        'axisMode', 'equal', ...
        'panelLayout', 'horizontal', ...
        'interpolationFactor', interpolationFactor, ...
        'mirrorCylindrical', true, ...
        'xlim', [-sampleRadius sampleRadius], ...
        'zlim', [0 sampleHeight], ...
        'boundaries', geometry.bnd);

    plotEnergies(plotting, obs, material, source.lambda);
    fprintf('Movies written to:\n  %s.gif\n  %s.avi\n',movieFile,movieFile);

end

% A check, that can be commented afterward!
cellVolumes = obs.dy .* (obs.dx(:) * obs.dz(:).');

integratedP = squeeze(sum(sum( ...
    PenergyDensity .* cellVolumes, 1), 2));

integratedS = squeeze(sum(sum( ...
    SenergyDensity .* cellVolumes, 1), 2));

integratedTotal = integratedP + integratedS;

figure;
plot(obs.t*1e6, integratedP, 'LineWidth', 1.5);
hold on;
plot(obs.t*1e6, integratedS, 'LineWidth', 1.5);
plot(obs.t*1e6, integratedTotal, 'k', 'LineWidth', 1.5);
yline(1, 'k--');
grid on;
xlabel('Time [microseconds]');
ylabel('Integrated energy');
legend('P', 'S', 'Total', 'Injected energy');
