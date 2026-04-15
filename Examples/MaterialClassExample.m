close all
clear variables
clc
%% define an acoustics material
for d = 2:3
    % geometry
    geometry = struct( 'dimension', d );

    SpectralLaw_def = {'exp','power_law','gaussian','triangular','low_pass'};
    for id = 1:5
        freq = 10; % in Hz
        matAcus = MaterialClass( geometry, ...
            freq, ...
            true, ...          % true for acoustics
            1, ...             % average wave velocity
            [0.1 0.2], ...     % coefficients of variation of kappa and rho.
            -0.5, ...          % correlation coefficient of kappa/rho
            SpectralLaw_def{id}, ...         % autocorrelation function
            0.1 );            % correlation length

        %% define an elastic material
        geometry = struct( 'dimension', d );
        freq = 10; % in Hz
        matElas = MaterialClass( geometry, ...
            freq, ...
            false, ...            % true for acoustics
            [6 6/sqrt(3)], ...    % defines the velocity of pressure waves and the shear waves
            [0.8 0.8 0.], ...     % defines the coefficients of variation of lambda, mu (Lamé coefficients) and rho (density), respectively.
            [0.1 0. 0.], ...      % defines the correlation coefficient between (lambda,mu), (lambda,rho), and (mu,rho), respectively
            SpectralLaw_def{id}, ...            % defines the autocorrelation function
            0.1);                 % defines the correlation length
        %% calc sigma and all the properties from the materials

        matAcus = matAcus.prepareSigma(matAcus,geometry.dimension);
        matElas = matElas.prepareSigma(matElas,geometry.dimension);

        disp(['Diffusion coefficient in Acoustic material: ', num2str(matAcus.Diffusivity)])
        disp(['Diffusion coefficient in Elastic material: ', num2str(matElas.Diffusivity)])
    end
end
%% plots

% sigma
a = figure;
matAcus.plotpolarsigma(a);
b = figure;
matElas.plotpolarsigma(b);

% PSDF
matAcus.PlotPSD;
matElas.PlotPSD;

% autocorrelation function
matAcus.PlotCorrelation;
matElas.PlotCorrelation;

% sigma not in polar
matAcus.plotsigma;
matElas.plotsigma;

%% create material with VonKarman ACF

matVK = MaterialClass();

matVK.d = 3;
matVK.acoustics = true;
matVK.v = 20;
matVK.rho = 10;
matVK.Frequency = 50;
matVK.coefficients_of_variation = [0.1 0.2];
matVK.correlation_coefficients =0;

matVK.SpectralLaw   = 'VonKarman';
matVK.SpectralParam = struct('nu',0.2);
matVK.CorrelationLength = matVK.v/matVK.Frequency;
matVK.getPSDF %it builds the power spectral density function
matVK.PlotCorrelation;
matVK.PlotPSD;
matVK = matVK.prepareSigma(matVK,matVK.d);

%% create material with monodisperse sphere/disk/rod ACF
% 3D - spheres
matMD3 = MaterialClass();
matMD3.d = 3;
matMD3.acoustics = true;
matMD3.v = 20;
matMD3.rho = 10;
matMD3.Frequency = 50;
matMD3.coefficients_of_variation = [0.1 0.2];
matMD3.correlation_coefficients = 0;
matMD3.SpectralLaw   = 'monodispersesphere';
matMD3.SpectralParam = struct('rhoS',0.2,'Diam',matMD3.v/matMD3.Frequency);
matMD3.getPSDF;   % calls MonoDisperseSphere internally (sets obj.R, obj.Phi, obj.CorrelationLength)
matMD3.PlotCorrelation;
matMD3.PlotPSD;
matMD3 = matMD3.prepareSigma(matMD3, matMD3.d);

% 2D - disks
matMD2 = MaterialClass();
matMD2.d = 2;
matMD2.acoustics = true;
matMD2.v = 20;
matMD2.rho = 10;
matMD2.Frequency = 50;
matMD2.coefficients_of_variation = [0.1 0.2];
matMD2.correlation_coefficients = 0;
matMD2.SpectralLaw   = 'monodispersesphere';
matMD2.SpectralParam = struct('rhoS',0.25,'Diam',matMD2.v/matMD2.Frequency);
matMD2.getPSDF;   % 2D PY OZ solver; also produces g(r)/S(k)/S2(r)/chi(r) figure
matMD2.PlotCorrelation;
matMD2.PlotPSD;
matMD2 = matMD2.prepareSigma(matMD2, matMD2.d);

% 1D - rods
matMD1 = MaterialClass();
matMD1.d = 1;
matMD1.acoustics = true;
matMD1.v = 20;
matMD1.Frequency = 50;
matMD1.coefficients_of_variation = [0.1 0.2];
matMD1.correlation_coefficients = 0;
matMD1.SpectralLaw   = 'monodispersesphere';
matMD1.SpectralParam = struct('rhoS',0.3,'Diam',matMD1.v/matMD1.Frequency);
matMD1.getPSDF;   % 1D PY OZ solver via HardBodyPY
matMD1.PlotCorrelation;
matMD1.PlotPSD;
return
%% pair correlation function and structure factor — standalone (HardBodyPY)
% MaterialClass.HardBodyPY returns g(r) and S(k) without constructing a full
% MaterialClass object. Useful for microstructural analysis or quick plotting.

phi_ex = 0.3;
D_ex   = matMD1.v/matMD1.Frequency;

[g1, S1, r1, k1] = MaterialClass.HardBodyPY(phi_ex, D_ex, 1);
[g2, S2, r2, k2] = MaterialClass.HardBodyPY(phi_ex, D_ex, 2);
[g3, S3, r3, k3] = MaterialClass.HardBodyPY(phi_ex, D_ex, 3);

figure('Name','HardBodyPY — g(r) comparison','Color','w','Position',[100 80 1100 380]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Hard-body pair correlation g(r),  \\phi = %.2f,  D = %.4g', phi_ex, D_ex), ...
    'FontSize',13,'FontWeight','bold');

nexttile;
plot(r1/D_ex, g1, 'b-', 'LineWidth', 1.8); hold on;
yline(1,'k--','LineWidth',0.8); xline(1,'r:','D','LineWidth',0.8);
xlabel('r / D'); ylabel('g(r)'); title('1D Hard Rods');
xlim([0 6]); grid on; box on;

nexttile;
plot(r2/D_ex, g2, 'b-', 'LineWidth', 1.8); hold on;
yline(1,'k--','LineWidth',0.8); xline(1,'r:','D','LineWidth',0.8);
xlabel('r / D'); ylabel('g(r)'); title('2D Hard Disks');
xlim([0 6]); grid on; box on;

nexttile;
plot(r3/D_ex, g3, 'b-', 'LineWidth', 1.8);
hold on;
yline(1,'k--','LineWidth',0.8); xline(1,'r:','D','LineWidth',0.8);
xlabel('r / D');
ylabel('g(r)');
title('3D Hard Spheres');
xlim([0 6]);
grid on;
box on;

figure('Name','HardBodyPY — S(k) comparison','Color','w','Position',[100 530 1100 380]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('Structure factor S(k),  \\phi = %.2f,  D = %.4g', phi_ex, D_ex), ...
    'FontSize',13,'FontWeight','bold');

nexttile;
plot(k1*D_ex, S1, 'r-', 'LineWidth', 1.8); hold on;
yline(1,'k--','LineWidth',0.8);
xlabel('k D'); ylabel('S(k)'); title('1D Hard Rods');
xlim([0 20]); grid on; box on;

nexttile;
plot(k2*D_ex, S2, 'r-', 'LineWidth', 1.8); hold on;
yline(1,'k--','LineWidth',0.8);
xlabel('k D'); ylabel('S(k)'); title('2D Hard Disks');
xlim([0 20]); grid on; box on;

nexttile;
% For d=3, k_grid = k_norm = k_phys*(D/2), so k_phys*D = 2*k_grid
plot(2*k3, S3, 'r-', 'LineWidth', 1.8); hold on;
yline(1,'k--','LineWidth',0.8);
xlabel('k D'); ylabel('S(k)'); title('3D Hard Spheres');
xlim([0 20]); grid on; box on;
%% create ''sphere'' sample

D = matMD1.v/matMD1.Frequency;
phi1 = 0.3;
L1 = [50 20 10]*2;
centersSphere = MaterialClass.CreateSphereComposite(L1,D,phi1);
MaterialClass.plot_map(centersSphere, D, L1);

phi2 = 0.2;
L2 = [100 50];
centersDisk = MaterialClass.CreateSphereComposite(L2,D,phi2);
MaterialClass.plot_map(centersDisk, D, L2);

phi3 = 0.35;
L3 = 100;
centersRod = MaterialClass.CreateSphereComposite(L3,D,phi3);
MaterialClass.plot_map(centersRod, D, L3);
% create a voxelization

resolution = 1;

% 3D — spheres  (L = [50 20 10], D = 1)
M3D = MaterialClass.VoxelizeDomain(centersSphere, D, L1, resolution);
fprintf('3D phi_vox = %.4f\n', mean(M3D(:)));

% 2D — disks  (L = [100 50], D = 1)
M2D = MaterialClass.VoxelizeDomain(centersDisk, D, L2, resolution);
fprintf('2D phi_vox = %.4f\n', mean(M2D(:)));

% 1D — rods  (L = 100, D = 1)
M1D = MaterialClass.VoxelizeDomain(centersRod, D, L3, resolution);
fprintf('1D phi_vox = %.4f\n', mean(M1D(:)));

% GetPSDFromImage: correlation function and PSDF from voxelized microstructures

% 3D — spheres
mat3D = MaterialClass();
mat3D.d = 3;
mat3D.GetPSDFromImage(M3D, resolution);

% 2D — disks
mat2D = MaterialClass();
mat2D.d = 2;
mat2D.GetPSDFromImage(M2D, resolution);

% 1D — rods
mat1D = MaterialClass();
mat1D.d = 1;
mat1D.GetPSDFromImage(M1D, resolution);

%% comparison: GetPSDFromImage (numerical) vs. MonoDisperseSphere PY theory

% Evaluate R(r/Lc) on a common normalised axis
r_cmp = linspace(0, 5, 500)';
R3_img = mat3D.R(r_cmp);
R3_py  = matMD3.R(r_cmp);
R2_img = mat2D.R(r_cmp);
R2_py  = matMD2.R(r_cmp);
R1_img = mat1D.R(r_cmp);
R1_py  = matMD1.R(r_cmp);

% Evaluate Phi(k*Lc) on a common normalised axis
k_cmp  = linspace(0, 6, 500)';
P3_img = mat3D.Phi(k_cmp);
P3_py  = matMD3.Phi(k_cmp);
P2_img = mat2D.Phi(k_cmp);
P2_py  = matMD2.Phi(k_cmp);
P1_img = mat1D.Phi(k_cmp);
P1_py  = matMD1.Phi(k_cmp);

% --- R(r/Lc) comparison ---
figure('Name','Image vs. PY — R(r/Lc)','Color','w','Position',[100 80 1200 380]);
tl = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
title(tl,'Normalised ACF: voxelised image (solid) vs. matMD PY theory (dashed)', ...
    'FontSize',13,'FontWeight','bold');

nexttile;
plot(r_cmp, R1_img, 'b-',  'LineWidth', 2.0);
hold on;
plot(r_cmp, R1_py,  'r--', 'LineWidth', 1.6);
yline(0,'k:','LineWidth',0.8);
xlabel('r / L_c');
ylabel('R(r)');
title(sprintf('1D Rods  (\\phi = %.2f)', phi3));
legend('Image','matMD (PY)','Location','northeast');
xlim([0 5]);
grid on;
box on;

nexttile;
plot(r_cmp, R2_img, 'b-',  'LineWidth', 2.0);
hold on;
plot(r_cmp, R2_py,  'r--', 'LineWidth', 1.6);
yline(0,'k:','LineWidth',0.8);
xlabel('r / L_c');
ylabel('R(r)');
title(sprintf('2D Disks  (\\phi = %.2f)', phi2));
legend('Image','matMD (PY)','Location','northeast');
xlim([0 5]);
grid on;
box on;

nexttile;
plot(r_cmp, R3_img, 'b-',  'LineWidth', 2.0);
hold on;
plot(r_cmp, R3_py,  'r--', 'LineWidth', 1.6);
yline(0,'k:','LineWidth',0.8);
xlabel('r / L_c');
ylabel('R(r)');
title(sprintf('3D Spheres  (\\phi = %.2f)', phi1));
legend('Image','matMD (PY)','Location','northeast');
xlim([0 5]);
grid on;
box on;

% --- Phi(k*Lc) comparison ---
figure('Name','Image vs. PY — Phi(k*Lc)','Color','w','Position',[100 510 1200 380]);
tl = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
title(tl,'PSDF: voxelised image (solid) vs. matMD PY theory (dashed)', ...
    'FontSize',13,'FontWeight','bold');

nexttile;
plot(k_cmp, P1_img, 'b-',  'LineWidth', 2.0); hold on;
plot(k_cmp, P1_py,  'r--', 'LineWidth', 1.6);
xlabel('k \cdot L_c');
ylabel('\Phi(k)');
title(sprintf('1D Rods  (\\phi = %.2f)', phi3));
legend('Image','matMD (PY)','Location','northeast');
xlim([0 6]);
grid on;
box on;

nexttile;
plot(k_cmp, P2_img, 'b-',  'LineWidth', 2.0);
hold on;
plot(k_cmp, P2_py,  'r--', 'LineWidth', 1.6);
xlabel('k \cdot L_c');
ylabel('\Phi(k)');
title(sprintf('2D Disks  (\\phi = %.2f)', phi2));
legend('Image','matMD (PY)','Location','northeast');
xlim([0 6]);
grid on;
box on;

nexttile;
plot(k_cmp, P3_img, 'b-',  'LineWidth', 2.0);
hold on;
plot(k_cmp, P3_py,  'r--', 'LineWidth', 1.6);
xlabel('k \cdot L_c');
ylabel('\Phi(k)');
title(sprintf('3D Spheres  (\\phi = %.2f)', phi1));
legend('Image','matMD (PY)','Location','northeast');
xlim([0 6]);
grid on;
box on;
