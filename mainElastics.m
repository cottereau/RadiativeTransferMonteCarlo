% This code works in 2D/3D. The initial condition is modeled as gaussian 
% with standard deviation source.lambda and uniformly-distributed angle(s).
% The initial direction is uniform.
clc

geometry = struct( 'dimension', 3 );

% Point source
source = struct( 'numberParticles', 1e6, ...
                 'polarization', 'S', ...
                 'lambda', 0.1 );

% observations
observation = struct('dr', 0.02, ...        % size of bins in space
                     'time', 0:0.2:4, ...  % observation times
                     'Ndir', 10 );         % number of bins for directions           

% material properties
% material.coefficients_of_variation defines the coefficients of variation
% of lambda, mu (Lamé coefficients) and rho (density), respectively.
% material.correlation_coefficients defines the correlation coefficient
% between (lambda,mu), (lambda,rho), and (mu,rho), respectively.

mat = MaterialClass(geometry.dimension);
mat.acoustics = false;
mat.vp = 6;
mat.vs = 6/sqrt(3);
mat.Frequency = 0.1;
mat.coefficients_of_variation = [0.8 0.8 0.];
mat.correlation_coefficients = [0.1 0. 0.];


% evaluate the Differential Scattering Cross-Sections
lc = 60;
mat.Exponential(lc);
mat.CalcSigma;
material = mat;
% % 2D
% material.sigma = {@(th) 1/2/pi*ones(size(th))*0.05*6, @(th) 1/2/pi*ones(size(th))*0.05*6; ...
%                   @(th) 1/2/pi*ones(size(th))*0.029*3.46, @(th) 1/2/pi*ones(size(th))*0.144*3.46};
% % 3D
% % material.sigma = {@(th) 1/8/pi*ones(size(th)), @(th) 2*sqrt(3)^3/20/pi*ones(size(th)); ...
% %                   @(th) 1/20/pi*ones(size(th)), @(th) 1/8/pi*ones(size(th))};

% radiative transfer solution - acoustic with boundaries
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% plotting output
plotting = struct( 'equipartition', true, ...
                   'movieTotalEnergy', false, ...
                   'movieDirectionalEnergy', false, ...
                   'timehistory', true, ...
                   'sensors', [1 0 0; 9 0 0; 9 0 1; 9 1 0]);

plotEnergies( plotting, obs, material, source.lambda );
