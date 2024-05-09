close all
clear all
clc

%% acoustic material
material = struct( 'acoustics', true, ...
                   'v', 1, ...
                   'Frequency', 10, ...                   
                   'spectralType', 'Ex', ...
                   'correlationLength', 10, ...
                   'coefficients_of_variation', [0.1 0.2], ...
                   'correlation_coefficients', -0.5 );
d = 3;

% evaluate the Differential Scattering Cross-Sections
dcs = DSCSClass(material,d);
dcs.CalcSigma;
dcs.Plotsigma
lc = dcs.CalcLc

% add Differential Scattering Cross-Sections to material
material.sigma = dcs.sigma;
