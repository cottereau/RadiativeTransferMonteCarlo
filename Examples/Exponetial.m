close all
clear all
clc

geometry = struct( 'type', 'fullspace', ...
                   'dimension', 3 );

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'position', [1 1 1]*0.1, ... 
                 'lambda', 0.01 );

% observations
observation = struct('dr', 0.01, ...        % size of bins in space
                     'time', 0:0.01:1, ...  % observation times
                     'Ndir', 10 );         % number of bins for directions           

% material properties
% material.coefficients_of_variation defines the coefficients of variation
% of kappa (bulk modulus) and rho (density), respectively.
% material.correlation_coefficients defines the correlation coefficient
% between kappa (bulk modulus) and rho (density).
mat = MaterialClass(geometry.dimension);
mat.acoustics = true;
mat.v = 1;
mat.Frequency = 10;
mat.coefficients_of_variation = [0.2 0.2];
mat.correlation_coefficients = -1;

% example monodisperse spheres
eta = 0.6;
D = 0.15;
%mat.MonoDisperseDisk(eta,D);

lc = 0.1;
mat.Exponential(lc);

mat.CalcSigma;
mat.PlotCorrelation;
mat.PlotPSD
%dcs.plotsigma
mat.plotpolarsigma
material = mat;
return
% a = figure(1);
% xlim([5 5.2])
% saveas(a,'ccoorr.png')
% 
% a = figure(2);
% xlim([0 60])
% saveas(a,'PSDFSphere.png')
% a = figure(3);
% title('Differential Scattering Cross-Section')
% saveas(a,'DSCSSphere.png')
% close all

% radiative transfer solution - acoustic with boundaries
switch geometry.type
    case 'fullspace'
        obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
    otherwise
        obs = radiativeTransferAcoustics( source, material, observation, geometry );
end

% plotting output
sensors = [1   1 1; 
           2   1 1;
           0.2 1 1;
           1.5 1 1;
           3.5 1 1]*0.1;

plotting = struct( 'equipartition', false, ...
                   'movieTotalEnergy', true, ...
                   'movieDirectionalEnergy', false, ...
                   'timehistory', true, ...
                   'sensors', sensors);

plotEnergies( plotting, obs, material, source.lambda )