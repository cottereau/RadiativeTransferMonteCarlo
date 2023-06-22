% This script contains some analytical validation cas studies for
% RadiativeTransferMonteCarlo code in fullspace  
%
% Our code is compared to the results described in the following papers:
% (1) J. C. J. Paasschens, Solution of the time-dependent Boltzmann equation, 
% Phys. Rev. E 56(1), pp. 1135-1141 (1997).
% (2) Hoshida (1991)
%
% Notes : 
% Other validation cases are expected to be added (higher dimensional and 
% anisotropic scattering in particular)
% 3D acoustic is currently under debug 

%% 2D Isotropic scattering (isotropic differential scattering cross-section)
% input data
titlecase = ['2D acoustic waves with isotropic differential scattering' ...
             ' cross section'];
disp(['(1) Testing ' titlecase ' ...']);
source      = struct( 'numberParticles', 1e6, ...
                      'position', [-10 0 0], ... 
                      'lambda', 0.001 );           
material    = struct( 'acoustics', true, ...
                      'v', 1, ...
                      'sigma', @(th) 1/4/pi*ones(size(th)));         
observation = struct( 'dr', 0.01, ...
                      'time', 0:0.1:15, ...
                      'Ndir', 100 );            
geometry    = struct( 'type', 'fullspace', ...
                      'size', [4 3 3], ...
                      'dimension', 2 );
inds = [40 80 120]; % index of the desired observation points
% running our code, Monte Carlo-based
obs = radiativeTransfer( source, material, observation, geometry );
Eus = obs.energyDensity.*obs.dr';
% computing Paasschens solution
EP = Paasschens_RTE_Unbounded( source, material, observation, geometry );
% running Hoshiba's Monte Carlo-based approach
EH = Hoshiba_RTE_Unbounded_MonteCarlo( obs.t, obs.r(inds), ...
                                          source, material, geometry, 10 );
EH = EH.*(2*pi*obs.dr(inds));
% visual comparison
figure; hold on; grid on; box on;
h1 = plot( obs.t, Eus(inds,:), '-k' );
h2 = plot( obs.t, EP(inds,:), '-b' );
h3 = plot( obs.t, EH, '-r' );
legend( [h1(1), h2(1), h3(1)], {'Monte Carlo (our code)','Analytical (Paasschens, 1997)', ...
        'Monte Carlo (Hoshiba 1991)'},'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);

%% 3D Isotropic scattering (isotropic differential scattering cross-section)
% input datatitlecase = ['3D acoustic waves with isotropic differential scattering' ...
             ' cross section'];
disp(['(2) Testing ' titlecase ' ...']);
source      = struct( 'numberParticles', 1e6, ...
                      'position', [-10 0 0], ... 
                      'lambda', 0.001 );           
material    = struct( 'acoustics', true, ...
                      'v', 1, ...
                      'sigma', @(th) 1/4/pi*ones(size(th)));         
observation = struct( 'dr', 0.05, ...
                      'time', 0:0.1:15, ...
                      'Ndir', 100 );            
geometry    = struct( 'type', 'fullspace', ...
                      'size', [4 3 3], ...
                      'dimension', 3 );
inds = [40 80 120]; % index of the desired observation points
% running our code, Monte Carlo-based
obs = radiativeTransfer( source, material, observation, geometry );
Eus = obs.energyDensity;
% computing Paasschens solution
EP = Paasschens_RTE_Unbounded( source, material, observation, geometry );
% % running Hoshiba's Monte Carlo-based approach
% EH = Hoshiba_RTE_Unbounded_MonteCarlo( obs.t, obs.r(inds), ...
%                                           source, material, geometry, 10 );
% visual comparison
figure; hold on; grid on; box on;
h1 = plot( obs.t, Eus(inds,:), '-k' );
h2 = plot( obs.t, EP(inds,:), '-b' );
%h3 = plot( obs.t, EH, '-r' );
legend( [h1(1), h2(1)], {'Monte Carlo (our code)','Analytical (Paasschens, 1997)'}, ...
        'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Energy density at different source-station distances')
title(titlecase);


%% 3D Anisotropic scattering (isotropic differential scattering cross-section)
% to be done
