close all
clearvars
clc

warning('To be done: 2D elastic differential scattering cross-sections to be investigated')

titlecase = '2D elastic case with anisotropic scattering';
disp(['Testing ' titlecase ' ...']);

geometry = struct( 'dimension', 2 );

source = struct('numberParticles', 5e6, ...
                'position', [0 0], ...
                'polarization', 'P', ...
                'lambda', 0.002);

% Random medium properties
freq = 2; % in Hz
material = MaterialClass(geometry, ...
    freq, ...            % frequency
    false, ...           % true for acoustics  
    [6 6/sqrt(3)], ...   % velocities of pressure and shear waves
    [0.1 0.1 0.05], ...  % coefficients of variation of lambda, mu (Lamé coefficients) and rho (density)
    [0. 0. 0.], ...         % correlation coefficients between (lambda,mu), (lambda,rho), and (mu,rho)
    'exp', ...           % autocorrelation function
    10);                 % correlation length

material = MaterialClass.prepareSigma( material, geometry.dimension );

% No intrinsic attenuation
material.Q = [Inf Inf];

observation = struct('x',0:1:10, ...
                     'y',[-pi pi], ...
                     'directions',[0 pi], ...
                     'time',0:0.01:1.5);

inds = [6 7 8]; % index of the desired observation points

% running our code, Monte Carlo-based
obs = radiativeTransfer( geometry, source, material, observation );

Ep = squeeze(obs.energyDensity(:,:,:,1));
Es = squeeze(obs.energyDensity(:,:,:,2));

% running Yoshimoto's Monte Carlo-based approach
EY = Comparison.randomWalkYoshimoto( geometry, source, material, observation, false );

% comparison of P energy densities
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Ep(inds,:), '-r');
h2 = semilogy( obs.t, EY(inds,:,1), '-b' );
set(gca, 'YScale', 'log');
legend( [h1(1), h2(1)], {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
xlabel('Lapse Time [s]')
ylabel('P-wave energy densities at different source-station distances')
title(titlecase);

% comparison of S energy densities
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Es(inds,:), '-r');
h2 = semilogy( obs.t, EY(inds,:,2), '-b' );
set(gca, 'YScale', 'log');
legend([h1(1), h2(1)],{'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
xlabel('Lapse Time [s]')
ylabel('S-wave energy densities at different source-station distances')
title(titlecase);