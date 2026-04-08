close all
clear variables
clc
% power_law 2D error

%% define an acoustics material
d = 3;
% geometry
geometry = struct( 'dimension', d );

SpectralLaw_def = {'exp','power_law','gaussian','triangular','low_pass'};
id = 5;
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
warning('Bug in Von Karman model')
if 0
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
    CorrelationLength       = matVK.v/matVK.Frequency;
    matVK.getPSDF %it builds the power spectral density function
    matVK.PlotCorrelation;
    matVK.PlotPSD;
    matVK = matVK.prepareSigma(matVK,matVK.d);
end
%% create material with  mono disperse sphere ACF
matMD3 = MaterialClass();

matMD3.d = 3;
matMD3.acoustics = true;
matMD3.v = 20;
matMD3.rho = 10;
matMD3.Frequency = 50;
matMD3.coefficients_of_variation = [0.1 0.2];
matMD3.correlation_coefficients =0;

matMD3.SpectralLaw   = 'monodispersesphere';
matMD3.SpectralParam = struct('rhoS',0.2,'Diam',matMD3.v/matMD3.Frequency);
matMD3.getPSDF %it builds the power spectral density function
matMD3.CalcLc;
matMD3.PlotCorrelation;
matMD3.PlotPSD;
matMD3 = matMD3.prepareSigma(matMD3,matMD3.d);

% 2D - disks
matMD2 = MaterialClass();

matMD2.d = 2;
matMD2.acoustics = true;
matMD2.v = 20;
matMD2.rho = 10;
matMD2.Frequency = 50;
matMD2.coefficients_of_variation = [0.1 0.2];
matMD2.correlation_coefficients =0;

matMD2.SpectralLaw   = 'monodispersesphere';
matMD2.SpectralParam = struct('rhoS',0.25,'Diam',matMD2.v/matMD2.Frequency);
matMD2.getPSDF %it builds the power spectral density function
matMD2.PlotCorrelation;
matMD2.PlotPSD;
matMD2 = matMD2.prepareSigma(matMD2,matMD2.d);
%% 'Imported'