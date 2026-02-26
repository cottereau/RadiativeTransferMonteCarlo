close all
clearvars
clc

titlecase = 'Comparison between different coordinate systems for 3D acoustic case with isotropic scattering';
disp(['Testing ' titlecase ' ...']);

%%%%%%%%%%%%%%%%%%%  Cylindrical %%%%%%%%%%%%%%%%%%%%%%

clear geometry material

geometry = struct( 'type', 'fullspace', ...
                   'dimension', 3 , ...
                   'frame', 'cylindrical' );
% Point source
source = struct( 'numberParticles', 1e6, ...
                 'type', 'point', ...         % 'point' (default) or 'plane'
                 'position', [0 0 0], ...    % always in cartesian frame
                 'direction', 'uniform',...  % with point source: 'uniform' (default) or 'outgoing'; with plane source: 1='x', 2='y', 3='z'
                 'lambda', 1e-4 );

observation = struct('x', 0:0.1:10, ...
                     'y', [-pi pi], ...
                     'z', [-1 1], ...
                     'directions', [0 pi], ...
                     'time', 0:0.05:20 );

material = MaterialClass.preset(1);
material.timeSteps = 0;

tic; obs_cyl = radiativeTransfer( geometry, source, material, observation ); toc;

%%%%%%%%%%%%%%%%%%%  Spherical %%%%%%%%%%%%%%%%%%%%%%
geometry.frame = 'spherical';

observation = struct('x', 0:0.1:10, ...
                     'y', [-pi pi], ...
                     'z', [-pi/2 pi/2], ...
                     'directions', [0 pi], ...
                     'time', observation.time );

tic; obs_sph = radiativeTransfer( geometry, source, material, observation ); toc;

%%%%%%%%%%%%%%%%%%%  Cartesian %%%%%%%%%%%%%%%%%%%%%%

geometry.frame = 'cartesian';

source.numberParticles = 2e6;

observation = struct('x', 0:0.1:10, ...
                     'y', [-1 1], ...
                     'z', [-1 1], ...
                     'directions', [0 pi], ...
                     'time', observation.time );

tic; obs_cart = radiativeTransfer( geometry, source, material, observation ); toc;

%%%%%%%%%%%%%%%%% Analytical %%%%%%%%%%%%%%%%%%%%%

EY = Comparison.randomWalkYoshimoto( geometry, source, material, observation );
[EP, Ediff] = Comparison.analyticalPaasschens( material, observation, geometry );

%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;

j = 50;

E_sph = squeeze(obs_sph.energyDensity(j,:,:));

E_cyl = squeeze(obs_cyl.energyDensity(j,:,:));

E_car = squeeze(obs_cart.energyDensity(j,:,:));

plot(obs_sph.t, E_sph);
plot(obs_sph.t, E_cyl);
plot(obs_sph.t, E_car);
plot(obs_sph.t, EY(j,:),'LineWidth',2);
plot(obs_sph.t, EP(j,:)); plot(obs_sph.t, Ediff(j,:));
set(gca,'YScale','log'); ylim([1e-4 2e-3]);
legend('spherical','cylindrical','cartesian','Yoshimoto','Paasschens','Diffusion','location','best');
xlabel('Time'); ylabel('Energy density'); grid on; box on;
