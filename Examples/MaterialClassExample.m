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
Dsphere = matMD3.v/matMD3.Frequency;
phi3 = 0.15;
matMD3.SpectralParam = struct('rhoS',phi3,'Diam',Dsphere);
matMD3.getPSDF;   % calls MonoDisperseSphere internally (sets obj.R, obj.Phi, obj.CorrelationLength)
matMD3.PlotCorrelation;
matMD3.PlotPSD;
matMD3 = matMD3.prepareSigma(matMD3, matMD3.d);

L1 = [50 50 50]*Dsphere;
centersSphere = MaterialClass.CreateSphereComposite(L1,Dsphere,phi3);
MaterialClass.plot_map(centersSphere, Dsphere, L1);

resolution = Dsphere/10;

M3D = MaterialClass.VoxelizeDomain(centersSphere, Dsphere, L1, resolution);
fprintf('3D phi_vox = %.4f\n', mean(M3D(:)));

% 3D — spheres
mat3D = MaterialClass();
mat3D.d = 3;
mat3D.GetPSDFromImage(M3D, resolution);

imageCorr =  mat3D.R( mat3D.r);
modelCorr = matMD3.R(matMD3.r);
%
figure
plot(mat3D.r,imageCorr)
hold on
plot(matMD3.r,modelCorr)
xlim([0 50*Dsphere])

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
phi2 = 0.15;
matMD2.SpectralParam = struct('rhoS',phi2,'Diam',Dsphere);
matMD2.getPSDF;   % 2D PY OZ solver; also produces g(r)/S(k)/S2(r)/chi(r) figure
matMD2.PlotCorrelation;
matMD2.PlotPSD;
matMD2 = matMD2.prepareSigma(matMD2, matMD2.d);


L2 = [50 50]*Dsphere;
centersSphere = MaterialClass.CreateSphereComposite(L2,Dsphere,phi2);
MaterialClass.plot_map(centersSphere, Dsphere, L2);

resolution = Dsphere/10;

M2D = MaterialClass.VoxelizeDomain(centersSphere, Dsphere, L2, resolution);
fprintf('2D phi_vox = %.4f\n', mean(M2D(:)));

% 2D — spheres
mat2D = MaterialClass();
mat2D.d = 2;
mat2D.GetPSDFromImage(M2D, resolution);

imageCorr =  mat2D.R( mat2D.r);
modelCorr = matMD2.R(matMD2.r);
%
figure
plot(mat2D.r,imageCorr)
hold on
plot(matMD2.r,modelCorr)
xlim([0 50*Dsphere])