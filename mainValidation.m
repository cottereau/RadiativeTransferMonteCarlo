% This script contains some analytical validation cas studies for
% RadiativeTransferMonteCarlo code in fullspace  

% Notes : 
% switch cases will be done when other cases are addes
% 3D acoustic is in debugging mode currently 
% Anisotropic scattering cases will be added

%% Isotropic scattering (differential scatteing cross-section independent of
% diection)

% Acoustic (scalar) waves in a 2D medium
% Point source
source = struct( 'numberParticles', 1e6, ...
                 'position', [-10 0 0], ... 
                 'lambda', 0.001 );
             
% material properties
material = struct( 'acoustics', true, ...
                   'v', 1, ...
                   'sigma', @(th) 1/4/pi*ones(size(th)) );
               
% observations
observation = struct('dr', 0.05, ...           % size of bins in space
                     'time', 0:0.1:15, ...     % observation times
                     'Ndir', 100, ...          % number of bins for directions
                     'sensors', [2 1.5 -1;
                                 2.9 1.5 -1]); % positions to plot directional energy            
 
geometry = struct( 'type', 'fullspace', ...
                   'size', [4 3 3], ...
                   'dimension', 2 );

% radiative transfer solution
obs = radiativeTransfer( source, material, observation, geometry );

% Analytical check using Paasschens 
EP = Paasschens_RTE_Unbounded(source, material, observation, geometry);

% Analytical check using Hoshiba's Monte Carlo-based approach (fixed obs point!)
EH = Hoshiba_RTE_Unbounded_MonteCarlo(obs.t,obs.r(80),source,material,geometry,30);

% Plot the result on a particular point
index = 80; % index of the desired observation point
figure; hold on;
plot(obs.t, obs.energyDensity(index,:)/(2*pi),'-k','LineWidth',2); % Why 2*pi is missing??!!
plot(obs.t, EP(index,:),'-b');
plot(obs.t,EH,'-r');
grid on; box on; legend('Monte Carlo (our code)', 'Paasschens', 'Monte Carlo (Hoshiba)');
xlabel('Lapse Time [s]');
ylabel(['Energy density at source-station distance of r = ', num2str(obs.r(index))])


%% Anisotropic scattering (to be done)