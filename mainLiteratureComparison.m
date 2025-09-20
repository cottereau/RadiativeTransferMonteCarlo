function mainLiteratureComparison(type)
% This function launches comparision cases for RadiativeTransferMonteCarlo
%
% syntax mainComparisonLiterature or mainComparisonLiterature(validationCase)
%
% when type='all' or empty, all validation cases are launched
%
% list of possible comparison cases with reference to literature
%  - '2dIsotropicAcoustic' (1,2,5)
%  - '3dIsotropicAcoustic' (1,2,5)
%  - '2dIsotropicElastic'  (3,5)
%  - '3dIsotropicElastic'  (4,5)
%  - '2dIsotropicAcoustic'
%  - '3dIsotropicAcoustic'
%  - '3dIsotropicElastic'
% Our code is compared to the results described in the following papers:
% (1) J. C. J. Paasschens. Solution of the time-dependent Boltzmann equation,
%     Phys. Rev. E 56(1), pp. 1135-1141 (1997).
% (2) M. Hoshiba. Simulation of multiple-scattered coda wave excitation
%     based on the energy conservation law. Phys. Earth Planet. Int. 67,
%     pp. 123-136 (1991).
% (3) H. Nakahara, K. Yoshimoto. Radiative transfer of elastic waves in
%     two-dimensional isotropic scattering media: semi-analytical approach
%     for isotropic source radiation. Earth Planets Space 63, pp. 459-468
%     (2011).
% (4) H. Sato. Multiple isotropic scattering model including P-S conversions
%     for the seismogram envelope formation. Geophys. J. Int 117,
%     pp. 487-494 (1994).
% (5) K. Yoshimoto, Monte Carlo simulation of seismogram envelopes in
%     scattering media. Journal of Geophysical Research: Solid Earth (2000)

if nargin==0
    mainLiteratureComparison('all')
else
    type = lower(type);
    switch type
        case 'all'
            mainLiteratureComparison('2dIsotropicAcoustic')
            mainLiteratureComparison('3dIsotropicAcoustic')
            mainLiteratureComparison('2dIsotropicElastic')
            mainLiteratureComparison('3disotropicElastic')
            mainLiteratureComparison('2dAnisotropicAcoustic')
            mainLiteratureComparison('3dAnisotropicAcoustic')
            %mainLiteratureComparison('2dAnisotropicElastic') % To be done
            mainLiteratureComparison('3dAnisotropicElastic')
            mainLiteratureComparison('3dIsotropicAcousticCoordinateSystems')
            mainLiteratureComparison('3dAnisotropicAcousticCoordinateSystems')
            
            %% 2D Isotropic scattering acoustic (isotropic differential scattering cross-section)
        case '2disotropicacoustic'
            titlecase = '2D acoustic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            
            % input data
            geometry = struct( 'dimension', 2 );
            
            source = struct( 'numberParticles', 1e6, ...
                             'position', [0 0], ...    
                             'lambda', 2e-4 );
            material = MaterialClass.preset(1);
            observation = struct('x', 0:0.03:9, ...
                'y', [-pi pi], ...
                'directions', [0 pi], ...
                'time', 0:0.05:20 );

            inds = [60 150 240]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry, source, material, observation );
            Eus = squeeze(obs.energyDensity);

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );

            % computing Paasschens solution
            [EP,Ediff] = Comparison.analyticalPaasschens( material, observation, geometry );

            % visual comparison
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Eus(inds,:), '-k' );
            h2 = semilogy( obs.t, EY(inds,:), '-r' );
            h3 = semilogy( obs.t, EP(inds,:), '-b' );
            h4 = semilogy( obs.t, Ediff(inds,:), ':r' );
            set(gca, 'YScale', 'log'); ylim([1e-5 1]);
            legend( [h1(1), h2(1), h3(1), h4(1)], ...
                {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
                'Analytical (Paasschens, 1997)', 'Diffusion approximation'}, ...
                'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('Integrated energy density at different source-station distances')
            title(titlecase);

            %% 3D Isotropic scattering (isotropic differential scattering cross-section)
        case '3disotropicacoustic'
            titlecase = '3D acoustic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            
            % input data
            geometry = struct( 'dimension', 3 );

            source = struct( 'numberParticles', 1e6, ...
                             'position', [0 0 0], ... 
                             'lambda', 2e-4 );

            material = MaterialClass.preset(1);

            observation = struct( 'x', 0:0.1:10, ...
                                  'y', [-pi pi], ...
                                  'z', [-pi/2 pi/2], ...
                                  'directions', [0 pi], ...
                                  'time', 0:0.05:20, ...   % observation times
                                  'Ndir', 10 );             % number of bins for directions

            inds = [20 50 80]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry, source, material, observation );
            Eus = squeeze(obs.energyDensity)/obs.dz;

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );

            % computing Paasschens solution
            [EP,Ediff] = Comparison.analyticalPaasschens( material, observation, geometry );

            % comparison
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Eus(inds,:), '-k' );
            h2 = semilogy( obs.t, EY(inds,:), '-r' );
            h3 = semilogy( obs.t, EP(inds,:), '-b' );
            h4 = semilogy( obs.t, Ediff(inds,:), ':r' );
            set(gca, 'YScale', 'log'); ylim([1e-5 1]);
            legend( [h1(1), h2(1), h3(1), h4(1)], ...
                {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
                'Analytical (Paasschens, 1997)', 'Diffusion approximation'}, ...
                'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('Integrated energy density at different source-station distances')
            title(titlecase);

            %% 2D Isotropic scattering elastic (Nakahara & Yoshimoto, 2011)
        case '2disotropicelastic'
            titlecase = '2D elastic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);

            geometry = struct( 'dimension', 2 );

            source = struct( 'numberParticles', 1e6, ...
                             'position', [0 0], ... 
                             'polarization', 'P', ...
                             'lambda', 0.002 );

            material = MaterialClass.preset(3);

            observation = struct('x', 0:0.1:20, ... % size of bins in space
                                 'y', [-pi pi], ...
                                 'directions', [0 pi], ...
                                 'time', 0:0.01:10 );

            d = geometry.dimension;
            vp = material.vp; vs = material.vs;
            K = vp/vs;
            Sigmapp = 0.2*vp; Sigmaps = 0.2*vp; Sigmap = Sigmapp + Sigmaps;
            Sigmasp = Sigmaps/((d-1)*K^d); Sigmass = Sigmap -Sigmasp;

            material.sigma = {@(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmapp, @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmaps; ...
                @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmasp, @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmass};

            inds = [20 50 80]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry, source, material, observation );

            Ep = squeeze(obs.energyDensity(:,:,:,1));
            Es = squeeze(obs.energyDensity(:,:,:,2));

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );
            %EtotY = sum(EY,3);

            % computing semi-analytical solution (Nakahara 2011, Sato 1994)
            % Eanalytical = Comparison.analyticalEnergyIsotropicElastic(geometry, material, observation, observation.r(inds), 0);

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

            %% 3D Isotropic scattering (isotropic differential scattering cross-section)
        case '3disotropicelastic'
            titlecase = '3D elastic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);

            geometry = struct( 'dimension', 3 );

            source = struct( 'numberParticles', 1e6, ...
                             'position', [0 0 0], ... 
                             'polarization', 'P', ...
                             'lambda', 0.002 );

            material = MaterialClass.preset(3);

            observation = struct('x', 0:0.1:20, ... % size of bins in space
                                 'y', [-pi pi], ...
                                 'z', [-pi/2 pi/2], ...
                                 'directions', [0 pi], ...
                                 'time', 0:0.01:10 );

            d = geometry.dimension;
            vp = material.vp; vs = material.vs;
            K = vp/vs;
            Sigmapp = 0.2*vp; Sigmaps = 0.2*vp; Sigmap = Sigmapp + Sigmaps;
            Sigmasp = Sigmaps/((d-1)*K^d); Sigmass = Sigmap -Sigmasp;

            material.sigma = {@(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmapp, @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmaps; ...
                @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmasp, @(th) 1/(2*pi^(d/2)/gamma(d/2))*ones(size(th))*Sigmass};

            Wsp = 0; % S to P wave energy ratio at the source

            inds = [20 50 80]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry, source, material, observation );

            Ep = squeeze(obs.energyDensity(:,:,:,1))/obs.dz;
            Es = squeeze(obs.energyDensity(:,:,:,2))/obs.dz;
            Etotus = Ep + Es;

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );
            EtotY = sum(EY,3);

            % computing semi-analytical solution (Nakahara 2011, Sato 1994)
            Eanalytical = Comparison.analyticalEnergyIsotropicElastic(geometry, material, observation, obs.x(inds), Wsp);

            % comparison of P energy densities
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Ep(inds,:), '-k' );
            h2 = semilogy( obs.t, EY(inds,:,1), '-b' );
            set(gca, 'YScale', 'log');
            legend( [h1(1), h2(1)], {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('P-wave energy densities at different source-station distances')
            title(titlecase);

            % comparison of S energy densities
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Es(inds,:), '-k' );
            h2 = semilogy( obs.t, EY(inds,:,2), '-b' );
            set(gca, 'YScale', 'log');
            legend( [h1(1), h2(1)], {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('S-wave energy densities at different source-station distances')
            title(titlecase);

            % comparison of total energy densities
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Etotus(inds,:), '-k' );
            h2 = semilogy( obs.t, EtotY(inds,:), '-r' );
            h3 = semilogy( obs.t, Eanalytical, '-b' );
            set(gca, 'YScale', 'log');
            legend( [h1(1), h2(1), h3(1)], {'Monte Carlo (our code)', ...
                'Monte Carlo (Yoshimoto 2000)','Analytical (Sato, 1994)'},'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('Energy densities at different source-station distances')
            currentYLim = ylim(gca);
            ylim(gca, [1e-5 currentYLim(2)]);
            title(titlecase);

            %% 2D Anisotropic scattering acoustic (anisotropic differential scattering cross-section)
        case '2danisotropicacoustic'
            titlecase = '2D acoustic case with anisotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            
            % input data
            geometry = struct( 'dimension', 2 );
            
            source = struct( 'numberParticles', 1e6, ...
                'position', [0 0], ...
                'lambda', 2e-4);

            material = MaterialClass.preset(1);
            % forward scattering regime
            material.sigma = {@(th) 1/4/pi*(1+4*cos(th).^4)};

            observation = struct( 'x', 0:0.03:9, ...
                                  'y', [-pi pi], ...
                                  'directions', [0 pi], ...
                                  'time', 0:0.02:20 );

            inds = [60 150 240]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry, source, material, observation );
            Eus = squeeze(obs.energyDensity);

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );

            % visual comparison
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Eus(inds,:), '-k' );
            h2 = semilogy( observation.time, EY(inds,:), '-r' );
            set(gca, 'YScale', 'log'); ylim([1e-5 1]);
            legend( [h1(1), h2(1)], ...
                {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('Integrated energy density at different source-station distances')
            title(titlecase);

            %% 3D Anisotropic scattering (anisotropic differential scattering cross-section)
        case '3danisotropicacoustic'
            titlecase = '3D acoustic case with anisotropic scattering';
            disp(['Testing ' titlecase ' ...']);
           
            % input data
            geometry = struct( 'dimension', 3 );

            source = struct( 'numberParticles', 1e6, ...
                'position', [0 0 0], ...
                'lambda', 2e-4 );

            %material = MaterialClass.preset(1);
            % forward scattering regime
            %material.sigma = {@(th) 1/4/pi*(1+4*cos(th).^4)};

            freq = 10; % in Hz
            material = MaterialClass( geometry, ...
                                      freq, ...
                                      true, ...          % true for acoustics
                                      1, ...             % average wave velocity
                                      [0.1 0.2], ...     % coefficients of variation of kappa and rho.
                                      -0.5, ...          % correlation coefficient of kappa/rho
                                      'exp', ...         % autocorrelation function
                                       0.1);              % correlation length

            observation = struct( 'x', 0:0.1:10, ...
                                  'y', [-pi pi], ...
                                  'z', [-pi/2 pi/2], ...
                                  'directions', [0 pi], ...
                                  'time', 0:0.02:20, ...  % observation times
                                  'Ndir', 10 );           % number of bins for directions

            inds = [20 50 80]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry, source, material, observation );
            Eus = squeeze(obs.energyDensity)/obs.dz;

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );

            % comparison
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Eus(inds,:), '-k' );
            h2 = semilogy( obs.t, EY(inds,:), '-r' );
            set(gca, 'YScale', 'log'); ylim([1e-5 1]);
            legend( [h1(1), h2(1)], ...
                {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('Integrated energy density at different source-station distances')
            title(titlecase);

            %% 2D Anisotropic scattering elastic
            % This case is to be further investigated since scattering
            % operators for 2D media are not given in Ryzhik et al, 1996
        case '2danisotropicelastic'
            titlecase = '2D elastic case with anisotropic scattering';
            disp(['Testing ' titlecase ' ...']);

            geometry = struct( 'dimension', 2 );

            source = struct( 'numberParticles', 1e6, ...
                             'position', [0 0], ...
                             'polarization', 'P', ...
                             'lambda', 0.002 );

            % The following setup generates a random medium favoring
            % forward scattering (large value for the normalized frequency)
            freq = 10; % in Hz
            material = MaterialClass( geometry, ...
                                      freq, ...
                                      false, ...            % true for acoustics
                                      [6 6/sqrt(3)], ...    % defines the velocity of pressure waves and the shear waves
                                      [0.8 0.8 0.], ...     % defines the coefficients of variation of lambda, mu (Lamé coefficients) and rho (density), respectively.
                                      [0.1 0. 0.], ...      % defines the correlation coefficient between (lambda,mu), (lambda,rho), and (mu,rho), respectively
                                      'exp', ...            % defines the autocorrelation function
                                       10);                 % defines the correlation length
            material = prepareSigma( material, geometry.dimension );

            observation = struct('x', 0:0.1:20, ... % size of bins in space
                                 'y', [-pi pi], ...
                                 'directions', [0 pi], ...
                                 'time', 0:0.01:10 );

            inds = [20 50 80]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry, source, material, observation );

            Ep = squeeze(obs.energyDensity(:,:,:,1));
            Es = squeeze(obs.energyDensity(:,:,:,2));

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation, false );

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

            %% 3D Anisotropic scattering (anisotropic differential scattering cross-section)
        case '3danisotropicelastic'
            titlecase = '3D elastic case with anisotropic scattering';
            disp(['Testing ' titlecase ' ...']);

            geometry = struct( 'dimension', 3 );

            source = struct( 'numberParticles', 1e6, ...
                             'position', [0 0 0], ...
                             'polarization', 'P', ...
                             'lambda', 0.002 );

            % The following setup favors a stochastic scattering regime
            freq = 10; % in Hz
            material = MaterialClass( geometry, ...
                                      freq, ...
                                      false, ...            % true for acoustics
                                      [6 6/sqrt(3)], ...    % defines the velocity of pressure waves and the shear waves
                                      [0.8 0.8 0.], ...     % defines the coefficients of variation of lambda, mu (Lamé coefficients) and rho (density), respectively.
                                      [0.1 0. 0.], ...      % defines the correlation coefficient between (lambda,mu), (lambda,rho), and (mu,rho), respectively
                                      'exp', ...            % defines the autocorrelation function
                                       0.1);                 % defines the correlation length
            material = prepareSigma( material, geometry.dimension );

            observation = struct('x', 0:0.1:20, ... % size of bins in space
                                 'y', [-pi pi], ...
                                 'z', [-pi/2 pi/2], ...
                                 'directions', [0 pi], ...
                                 'time', 0:0.01:10 );

            inds = [20 50 80]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry, source, material, observation );

            Ep = squeeze(obs.energyDensity(:,:,:,1))/obs.dz;
            Es = squeeze(obs.energyDensity(:,:,:,2))/obs.dz;
            Etotus = Ep + Es;

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation, false );
            EtotY = sum(EY,3);

            % comparison of P energy densities
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Ep(inds,:), '-k' );
            h2 = semilogy( obs.t, EY(inds,:,1), '-b' );
            set(gca, 'YScale', 'log');
            legend( [h1(1), h2(1)], {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('P-wave energy densities at different source-station distances')
            title(titlecase);

            % comparison of S energy densities
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Es(inds,:), '-k' );
            h2 = semilogy( obs.t, EY(inds,:,2), '-b' );
            set(gca, 'YScale', 'log');
            legend( [h1(1), h2(1)], {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('S-wave energy densities at different source-station distances')
            title(titlecase);

            % comparison of total energy densities
            figure; hold on; grid on; box on;
            h1 = semilogy( obs.t, Etotus(inds,:), '-k' );
            h2 = semilogy( obs.t, EtotY(inds,:), '-r' );
            set(gca, 'YScale', 'log');
            legend( [h1(1), h2(1)], {'Monte Carlo (our code)','Monte Carlo (Yoshimoto 2000)'},'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('Energy densities at different source-station distances')
            currentYLim = ylim(gca);
            ylim(gca, [1e-5 currentYLim(2)]);
            title(titlecase);

        %% 3D Isotropic scattering acoustic (comparison between cylindrical, spherical, cartesian versions of the code)
        case '3disotropicacousticcoordinatesystems'
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
            
            tic; obs_cyl = radiativeTransferUnbounded( geometry, source, material, observation ); toc;
            
            %%%%%%%%%%%%%%%%%%%  Spherical %%%%%%%%%%%%%%%%%%%%%%
            geometry.frame = 'spherical';
            
            observation = struct('x', 0:0.1:10, ...
                                 'y', [-pi pi], ...
                                 'z', [-pi/2 pi/2], ...
                                 'directions', [0 pi], ...
                                 'time', observation.time );
            
            tic; obs_sph = radiativeTransferUnbounded( geometry, source, material, observation ); toc;
            
            %%%%%%%%%%%%%%%%%%%  Cartesian %%%%%%%%%%%%%%%%%%%%%%
            
            geometry.frame = 'cartesian';
            
            source.numberParticles = 2e6;
            
            observation = struct('x', 0:0.1:10, ...
                                 'y', [-1 1], ...
                                 'z', [-1 1], ...
                                 'directions', [0 pi], ...
                                 'time', observation.time );
            
            tic; obs_cart = radiativeTransferUnbounded( geometry, source, material, observation ); toc;
            
            %%%%%%%%%%%%%%%%% Analytical %%%%%%%%%%%%%%%%%%%%%
            
            EY = Comparison.randomWalkYoshimoto_beta( geometry, source, material, observation );
            [EP,Ediff] = Comparison.analyticalPaasschens( material, observation, geometry );
            
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
            plot(obs_sph.t, EP(j,:)); plot(obs_sph.t, Ediff(j,:));
            set(gca,'YScale','log'); ylim([1e-4 2e-3]);
            legend('spherical','cylindrical','cartesian','Yoshimoto','Paasschens','Diffusion','location','best');
            xlabel('Time'); ylabel('Energy density'); grid on; box on;

        %% 3D Anisotropic scattering acoustic (comparison between cylindrical, spherical, cartesian versions of the code)
        case '3danisotropicacousticcoordinatesystems'
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
            
            tic; obs_cyl = radiativeTransferUnbounded( geometry, source, material, observation ); toc;
            
            %%%%%%%%%%%%%%%%%%%  Spherical %%%%%%%%%%%%%%%%%%%%%%
            
            geometry.frame = 'spherical';
            
            observation = struct('x', 0:0.1:10, ...
                                 'y', [-pi pi], ...
                                 'z', [-pi/2 pi/2], ...
                                 'directions', [0 pi], ...
                                 'time', observation.time );
            
            tic; obs_sph = radiativeTransferUnbounded( geometry, source, material, observation ); toc;
            
            %%%%%%%%%%%%%%%%%%%  Cartesian %%%%%%%%%%%%%%%%%%%%%%
            
            geometry.frame = 'cartesian';
            
            source.numberParticles = 4e6;
            
            observation = struct('x', 0:0.1:10, ...
                                 'y', [-1 1], ...
                                 'z', [-1 1], ...
                                 'directions', [0 pi], ...
                                 'time', observation.time );
            
            tic; obs_cart = radiativeTransferUnbounded( geometry, source, material, observation ); toc;
            
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

        otherwise
            error('unknown validation case')
    end
end