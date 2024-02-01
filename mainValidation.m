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

%% 2D Isotropic scattering acoustic (isotropic differential scattering cross-section)
% input data
titlecase = ['2D acoustic waves with isotropic differential scattering' ...
             ' cross section'];
disp(['(1) Testing ' titlecase ' ...']);
source      = struct( 'numberParticles', 1e6, ...
                      'position', [-10 0 0], ... 
                      'lambda', 0.001 );           
material    = struct( 'acoustics', true, ...
                      'v', 2, ...
                      'sigma', @(th) 1/4/pi*ones(size(th)));         
observation = struct( 'dr', 0.01, ...
                      'time', 0:0.1:15, ...
                      'Ndir', 100 );            
geometry    = struct( 'type', 'fullspace', ...
                      'size', [4 3 3], ...
                      'dimension', 2 );
inds = [40 80 120]; % index of the desired observation points
% running our code, Monte Carlo-based
%obs = radiativeTransferAcoustics( source, material, observation, geometry );
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
Eus = (obs.Ei+obs.Ec).*obs.dr';
% computing Paasschens solution
[EP,Ediff] = Paasschens_RTE_Unbounded( source, material, observation, geometry );
% running Hoshiba's Monte Carlo-based approach
EH = Hoshiba_RTE_Unbounded_MonteCarlo( obs.t, obs.r(inds), ...
                                          source, material, geometry, 10 );
EH = EH.*(2*pi*obs.dr(inds));
% visual comparison
figure; hold on; grid on; box on;
h1 = plot( obs.t, Eus(inds,:), '-k' );
h2 = plot( obs.t, EP(inds,:), '-b' );
h3 = plot( obs.t, EH, '-r' );
h4 = plot( obs.t, Ediff(inds,:), ':r' );
legend( [h1(1), h2(1), h3(1), h4(1)], ...
 {'Monte Carlo (our code)','Analytical (Paasschens, 1997)', ...
  'Monte Carlo (Hoshiba 1991)', 'diffusion approximation'}, ...
'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);

%% 3D Isotropic scattering (isotropic differential scattering cross-section)
% input data
titlecase = ['3D acoustic waves with isotropic differential scattering' ...
             ' cross section'];
disp(['(2) Testing ' titlecase ' ...']);
source      = struct( 'numberParticles', 1e6, ...
                      'position', [-10 0 0], ... 
                      'lambda', 0.001 );           
material    = struct( 'acoustics', true, ...
                      'v', 2, ...
                      'sigma', @(th) 1/10/pi*ones(size(th)));         
observation = struct( 'dr', 0.05, ...
                      'time', 0:0.1:30, ...
                      'Ndir', 100 );            
geometry    = struct( 'type', 'fullspace', ...
                      'size', [4 3 3], ...
                      'dimension', 3 );
inds = [40 80 120]; % index of the desired observation points
% running our code, Monte Carlo-based
% obs = radiativeTransferAcoustics( source, material, observation, geometry );
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
Eus = (obs.Ei+obs.Ec).*obs.dr';
% computing Paasschens solution
[EP,Ediff] = Paasschens_RTE_Unbounded( source, material, observation, geometry );
% % running Hoshiba's Monte Carlo-based approach
% EH = Hoshiba_RTE_Unbounded_MonteCarlo( obs.t, obs.r(inds), ...
%                                           source, material, geometry, 10 );
% visual comparison
figure; hold on; grid on; box on;
h1 = plot( obs.t, Eus(inds,:), '-k' );
h2 = plot( obs.t, EP(inds,:), '-b' );
%h3 = plot( obs.t, EH, '-r' );
h4 = plot( obs.t, Ediff(inds,:), ':r' );
legend( [h1(1), h2(1), h4(1)], ...
 {'Monte Carlo (our code)','Analytical (Paasschens, 1997)', ...
  'diffusion approximation'}, ...
'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);

%% 2D Isotropic scattering elastic (Nakahara & Yoshimoto, 2011)
titlecase = ['2D elastic waves with isotropic differential scattering' ...
             ' cross sections'];
% Load the results of the paper
load('Nakahara_Yoshimoto_2011_Fig4_left.mat')

% Run our code
geometry = struct( 'dimension', 2 );

source = struct( 'numberParticles', 1e6, ...
                 'polarization', 'S', ...
                 'lambda', 0.001 );

observation = struct('dr', 0.4, ...        % size of bins in space
                     'time', 0:0.0167:34.13, ...  % observation times
                     'Ndir', 10 );         % number of bins for directions           
 
% material properties
% material.coefficients_of_variation defines the coefficients of variaiton
% of lambda, mu (Lamé coefficients) and rho (density), respectively.
% material.correlation_coefficients defines the correlation coefficient
% between (lambda,mu), (lambda,rho), and (mu,rho), respectively.
material = struct( 'acoustics', false, ...
                   'vp', 6, ...
                   'vs', 3.46);
% Values taken from (Nakahara & Yoshimoto, 2011, Section 3.1)
K = material.vp/material.vs;
gpp = 0.05; gps = 0.05; gp = gpp + gps;
gsp = gps/K; gss = (material.vp*gp - material.vs*gsp)/material.vs;
material.sigma = {@(th) 1/2/pi*ones(size(th))*gpp*material.vp, @(th) 1/2/pi*ones(size(th))*gps*material.vp; ...
                  @(th) 1/2/pi*ones(size(th))*gsp*material.vs, @(th) 1/2/pi*ones(size(th))*gss*material.vs};

obsS = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

% Epc = obs.Ec(:,:,1); Epi = obs.Ei(:,:,1); Ep = Epc+Epi; Esc = obs.Ec(:,:,2); Esi = obs.Ei(:,:,2); Es = Esc+Esi;
EpcP = obsP.Ec(:,:,1); EpiP = obsP.Ei(:,:,1); EpP = EpcP+EpiP; EscP = obsP.Ec(:,:,2); EsiP = obsP.Ei(:,:,2); EsP = EscP+EsiP;
EpcS = obsS.Ec(:,:,1); EpiS = obsS.Ei(:,:,1); EpS = EpcS+EpiS; EscS = obsS.Ec(:,:,2); EsiS = obsS.Ei(:,:,2); EsS = EscS+EsiS;
Ep = (1*EpP+1.5*K^5*EpS)/(1+1.5*K^5); Es = (1*EsP+1.5*K^5*EsS)/(1+1.5*K^5);
% Normalized P-wave energy in terms of normalized time
figure; semilogy(obs.t*gp*material.vp,2*pi*Ep(26,:)/(1/(1+1.5*K^5)*(gp*material.vp)^2),'-b','linewidth',2);
xlabel('Normalized Time [-]'); ylabel('Normalized P-wave energy [-]');
xlim([0 6]); ylim([0.001 100]); hold on; plot(energy_P(:,1),energy_P(:,2),'-r','linewidth',2);
yticks([0.001 0.01 0.1 1 10 100]); yticklabels({'0.001','0.01','0.1','1','10','100'});

% Normalized P-wave energy in terms of normalized time
figure; semilogy(obs.t*gp*material.vp,2*pi*Es(26,:)/(1/(1+1.5*K^5)*(gp*material.vp)^2),'-b','linewidth',2); 
xlabel('Normalized Time [-]'); ylabel('Normalized S-wave energy [-]'); 
xlim([0 6]); ylim([0.001 100]); hold on; plot(energy_S(:,1),energy_S(:,2),'-r','linewidth',2);
yticks([0.001 0.01 0.1 1 10 100]); yticklabels({'0.001','0.01','0.1','1','10','100'});

%% 3D Anisotropic scattering (isotropic differential scattering cross-section)
% to be done



return

%% Generic Random Input Parameters
% input data
titlecase = 'Generic Random Inputs ...';
disp(['(3) Testing ' titlecase ' ...']);
source      = struct( 'numberParticles', 1e6, ...
                      'position', [-10 0 0], ... 
                      'lambda', 0.001 );                 
randpars    = struct( 'normfreq', 1, ...
                      'corr_distance', 1, ...
                      'PSDF_type','exp', ...
                      'corr_matrix',[0.1 0; 0 0.1]);
if ~issymmetric(randpars.corr_matrix)
    disp('corr_matirx field should be a symmetric matrix !')
end
% sigma obtained in the following line is normalized by w (angular freq)
% for sigmaSP and sigmaSS, multiplication by the identity matrix is not done 
% everywhere in this version
sigma = PSDF2sigma(acoustics, material, randpars);

material    = struct( 'acoustics', true, ...
                      'v', 1, ...
                      'sigma', sigma);  

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
[EP,Ediff] = Paasschens_RTE_Unbounded( source, material, observation, geometry );
% running Hoshiba's Monte Carlo-based approach
EH = Hoshiba_RTE_Unbounded_MonteCarlo( obs.t, obs.r(inds), ...
                                          source, material, geometry, 10 );
EH = EH.*(2*pi*obs.dr(inds));
% visual comparison
figure; hold on; grid on; box on;
h1 = plot( obs.t, Eus(inds,:), '-k' );
h2 = plot( obs.t, EP(inds,:), '-b' );
h3 = plot( obs.t, EH, '-r' );
h4 = plot( obs.t, Ediff(inds,:), ':r' );
legend( [h1(1), h2(1), h3(1), h4(1)], ...
 {'Monte Carlo (our code)','Analytical (Paasschens, 1997)', ...
  'Monte Carlo (Hoshiba 1991)', 'diffusion approximation'}, ...
'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);

