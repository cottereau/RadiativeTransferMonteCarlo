close all
clear all
clc

titlecase = 'Comparison between different coordinate systems for 3D acoustic case with anisotropic scattering';
disp(['Testing ' titlecase ' ...']);

%%%%%%%%%%%%%%%%%%%  Cylindrical %%%%%%%%%%%%%%%%%%%%%%
clear geometry material

geometry = struct( 'type', 'fullspace', ...
    'dimension', 3 , ...
    'frame', 'cylindrical' );

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
cv = 10; ell = 0.5;
freq = 1;
velo = 01;
material = MaterialClass( geometry, ...
    freq, ...
    true, ...          % true for acoustics
    velo, ...             % average wave velocity
    [cv cv]/100, ...     % coefficients of variation of kappa and rho.
    0, ...          % correlation coefficient of kappa/rho
    'Gauss', ...         % autocorrelation function
    ell);              % correlation length
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

source.numberParticles = 4e6;

observation = struct('x', 0:0.1:10, ...
    'y', [-1 1], ...
    'z', [-1 1], ...
    'directions', [0 pi], ...
    'time', observation.time );

tic; obs_cart = radiativeTransfer( geometry, source, material, observation ); toc;

%%%%%%%%%%%%%%%%% Analytical %%%%%%%%%%%%%%%%%%%%%

EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );

%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%
figure; j = 50;

E_sph = squeeze(obs_sph.energyDensity(j,:,:)) ./ obs_sph.dz;

% exact cylindrical correction using bin edges:
rL = obs_cyl.binX(j); rR = obs_cyl.binX(j+1);
int_r2 = obs_cyl.dx(j);                % what the code uses for normalization
int_r  = 0.5*(rR^2 - rL^2);            % what cylindrical volume needs
E_cyl = squeeze(obs_cyl.energyDensity(j,:,:)) .* ( int_r2 / (int_r * obs_cyl.dz) );

E_car = squeeze(obs_cart.energyDensity(j,:,:)) ./ obs_cart.dz;

plot(obs_sph.t, E_sph); hold on
plot(obs_sph.t, E_cyl);
plot(obs_sph.t, E_car);
plot(obs_sph.t, EY(j,:),'LineWidth',2);
set(gca,'YScale','log');
legend('spherical','cylindrical','cartesian','Yoshimoto','location','best');
xlabel('Time'); ylabel('Energy density'); grid on; box on;