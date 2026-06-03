%% create material with monodisperse sphere/disk/rod ACF
ccc

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
Dsphere = matMD3.v/matMD3.Frequency;
phi3 = 0.2;
matMD3.SpectralParam = struct('rhoS',phi3,'Diam',Dsphere);
matMD3.getPSDF;   
matMD3 = matMD3.prepareSigma(matMD3, matMD3.d);


L1 = [30 30 30]*Dsphere;
centersSphere = MaterialClass.CreateSphereComposite(L1,Dsphere,phi3);
MaterialClass.plot_map(centersSphere, Dsphere, L1);

resolution = Dsphere/5;

M3D = MaterialClass.VoxelizeDomain(centersSphere, Dsphere, L1, resolution);
fprintf('3D phi_vox = %.4f\n', mean(M3D(:)));
MaterialClass.PlotVoxel(M3D, resolution);
ylim([0 5])
% 3D — spheres
clear mat3D
mat3D = MaterialClass();
mat3D.d = 3;
mat3D.GetPSDFromImage(M3D, resolution);

% Compare in physical r/D units (independent of Lc normalisation)
figure;
r_phys = linspace(0, 6*Dsphere, 500);
plot(r_phys/Dsphere, matMD3.R(r_phys/matMD3.CorrelationLength), 'b-', 'LineWidth', 2)
hold on
plot(r_phys/Dsphere, mat3D.R(r_phys/mat3D.CorrelationLength), 'r--', 'LineWidth', 2)
xlabel('r / D'); ylabel('R(r)')
legend('Analytical (PY)', 'Microstructure')
title('3D — ACF comparison (r/D)')
xlim([0 6]); grid on; box on

% Compare physical PSD in k*D units: Phi_phys(k) = Lc^3 * Phi_norm(k*Lc)
figure;
k_D = linspace(1e-3, 20, 500);          % k*D, dimensionless
k_phys = k_D / Dsphere;                  % physical k [1/D units]
Lc_MD3 = matMD3.CorrelationLength;
Lc_3D  = mat3D.CorrelationLength;
PSD_MD3 = (Lc_MD3/Dsphere)^3 * matMD3.Phi(k_phys * Lc_MD3);
PSD_3D  = (Lc_3D/Dsphere)^3  * mat3D.Phi( k_phys * Lc_3D);
plot(k_D, PSD_MD3, 'b-',  'LineWidth', 2)
hold on
plot(k_D, PSD_3D,  'r--', 'LineWidth', 2)
xlabel('k \cdot D'); ylabel('\Phi(k) / D^3')
legend('Analytical (PY)', 'Microstructure')
title('3D — PSD comparison (k\cdotD)')
xlim([0 20]); grid on; box on
%%

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
Dsphere = matMD2.v/matMD2.Frequency;
phi2 = 0.2;
matMD2.SpectralParam = struct('rhoS',phi2,'Diam',Dsphere);
matMD2.getPSDF;
matMD2 = matMD2.prepareSigma(matMD2, matMD2.d);

L2 = 2*[50 50]*Dsphere;
centersSphere = MaterialClass.CreateSphereComposite(L2,Dsphere,phi2);
MaterialClass.plot_map(centersSphere, Dsphere, L2);

resolution = Dsphere/5;

M2D = MaterialClass.VoxelizeDomain(centersSphere, Dsphere, L2, resolution);
fprintf('2D phi_vox = %.4f\n', mean(M2D(:)));
MaterialClass.PlotVoxel(M2D, resolution);

% 2D — spheres
clear mat2D
mat2D = MaterialClass();
mat2D.d = 2;
mat2D.GetPSDFromImage(M2D, resolution);

% Compare in physical r/D units
figure;
r_phys2 = linspace(0, 6*Dsphere, 500);
plot(r_phys2/Dsphere, matMD2.R(r_phys2/matMD2.CorrelationLength), 'b-', 'LineWidth', 2)
hold on
plot(r_phys2/Dsphere, mat2D.R(r_phys2/mat2D.CorrelationLength), 'r--', 'LineWidth', 2)
xlabel('r / D'); ylabel('R(r)')
legend('Analytical (PY)', 'Microstructure')
title('2D — ACF comparison (r/D)')
xlim([0 6]); grid on; box on

% Compare physical PSD in k*D units: Phi_phys(k) = Lc^2 * Phi_norm(k*Lc)  (d=2)
figure;
k_D2 = linspace(1e-3, 20, 500);
k_phys2 = k_D2 / Dsphere;
Lc_MD2 = matMD2.CorrelationLength;
Lc_2D  = mat2D.CorrelationLength;
PSD_MD2 = (Lc_MD2/Dsphere)^2 * matMD2.Phi(k_phys2 * Lc_MD2);
PSD_2D  = (Lc_2D/Dsphere)^2  * mat2D.Phi( k_phys2 * Lc_2D);
plot(k_D2, PSD_MD2, 'b-',  'LineWidth', 2)
hold on
plot(k_D2, PSD_2D,  'r--', 'LineWidth', 2)
xlabel('k \cdot D'); ylabel('\Phi(k) / D^2')
legend('Analytical (PY)', 'Microstructure')
title('2D — PSD comparison (k\cdotD)')
xlim([0 20]); grid on; box on