close all
clear variables
clc

geometry = struct( 'type', 'fullspace', ...
                   'dimension', 3 );

% Point source
source = struct( 'numberParticles', 10e6, ...
                 'position', [0 0 0], ... 
                 'lambda', 0.001 );

% observations
observation = struct('dr', 0.0025, ...        % size of bins in space
                     'time', 0:0.025:2, ...  % observation times
                     'Ndir', 20 );         % number of bins for directions           

% material properties
% material.coefficients_of_variation defines the coefficients of variation
% of kappa (bulk modulus) and rho (density), respectively.
% material.correlation_coefficients defines the correlation coefficient
% between kappa (bulk modulus) and rho (density).
mat = MaterialClass(geometry.dimension);
mat.acoustics = true;
mat.v = 1;
mat.Frequency = 1;
mat.coefficients_of_variation = [0.1 0.1];
mat.correlation_coefficients = 1;

% example of an image
Im=imread('be1.png');
dx = 0.1e-3;
imshow(Im)
dir = 1;

%dx = mean(diff(x));
mat.GetPSDFromImage(Im,dir,dx,dx);
return
mat.CalcSigma;
mat.PlotPSD

%dcs.plotsigma
mat.plotpolarsigma
material = mat;
% 
a = figure(1);
saveas(a,'Sample.png')
a = figure(2);
xlim([0 10])
saveas(a,'PSDF.png')
a = figure(3);
title('Differential Scattering Cross-Section')
saveas(a,'DSCS.png')

% radiative transfer solution - acoustic with boundaries
switch geometry.type
    case 'fullspace'
        obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
    otherwise
        obs = radiativeTransferAcoustics( source, material, observation, geometry );
end

% plotting output
sensors = [0 0 0;
    0 0 40;
    0 0 80;
    0 0 120;
    0 0 160;
    0 0 200;
    0 0 240;
    0 0 280;
    0 0 320;
    0 0 360;
    0 0 400]*1e-3;

plotting = struct( 'equipartition', false, ...
                   'movieTotalEnergy', true, ...
                   'movieDirectionalEnergy', false, ...
                   'timehistory', true, ...
                   'sensors', sensors);

plotEnergies( plotting, obs, material, source.lambda )