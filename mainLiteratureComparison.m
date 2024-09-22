function mainLiteratureComparison(type)
% This function launches comparision cases for RadiativeTransferMonteCarlo
%
% syntax mainComparisonLiterature or mainComparisonLiterature(validationCase)
%
% when type='all' or empty, all validation cases are launched
%
% list of possible comparision cases with reference to literature
%  - '2dIsotropicAcoustic' (1,2,5)
%  - '3dIsotropicAcoustic' (1,2,5)
%  - '2dIsotropicElastic'  (3,5)
%  - '3dIsotropicElastic'  (4,5)
%
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

% with no argument, launch all validation cases


% Bug reports: 
% Analytical comparison in 2D does not work at the moment (corresponding
% line is commented)

if nargin==0
    mainLiteratureComparison('all')
else
    type = lower(type);
    switch type
        case 'all'
            mainLiteratureComparison('2dIsotropicAcoustic')
            mainLiteratureComparison('3dIsotropicAcoustic')
            mainLiteratureComparison('2dIsotropicElastic')
            mainLiteratureComparison('3disotropicelastic')

        %% 2D Isotropic scattering acoustic (isotropic differential scattering cross-section)
        case '2disotropicacoustic'
            titlecase = '2D acoustic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            % input data
            geometry = struct( 'dimension', 2 );
            source = struct( 'numberParticles', 1e6, ...
                             'lambda', 2e-4 );
            material = MaterialClass.preset(1);
            observation = struct('r', 0:0.03:9, ...
                                 'azimuth', [-pi pi], ... 
                                 'directions', [0 pi], ...             
                                 'time', 0:0.05:20 );
            
            inds = [60 150 240]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
            Eus = squeeze(obs.energyDensity);

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto( geometry, source, material, observation, false );

            % computing Paasschens solution
            [EP,Ediff] = Comparison.analyticalPaasschens( material, observation, geometry );

            % visual comparison
            figure; hold on; grid on; box on;
            h1 = plot( obs.t, Eus(inds,:), '-k' );
            h2 = plot( obs.t, EY(inds,:), '-r' );
            h3 = plot( obs.t, EP(inds,:), '-b' );
            h4 = plot( obs.t, Ediff(inds,:), ':r' );
            legend( [h1(1), h2(1), h3(1), h4(1)], ...
                {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
                 'Analytical (Paasschens, 1997)', 'diffusion approximation'}, ...
                 'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('Integrated energy density at different source-station distances')
            title(titlecase);

        %% 3D Isotropic scattering (isotropic differential scattering cross-section)
        case '3disotropicacoustic'
            titlecase = '3D acoustic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            % input data
            source = struct( 'numberParticles', 1e6, ...
                             'lambda', 2e-4 );
            
            material = MaterialClass.preset(1);

            observation = struct( 'r', 0:0.1:10, ...
                                  'azimuth', [-pi pi], ... 
                                  'elevation', [-pi/2 pi/2], ... 
                                  'directions', [0 pi], ...             
                                  'time', 0:0.05:20, ...   % observation times
                                  'Ndir', 10 );             % number of bins for directions

            geometry = struct( 'type', 'fullspace', ...
                                'dimension', 3 );
            
            inds = [20 50 80]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
            Eus = squeeze(obs.energyDensity)/obs.dphi;

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto( geometry, source, material, observation, false );

            % computing Paasschens solution
            [EP,Ediff] = Comparison.analyticalPaasschens( material, observation, geometry );

            % comparison
            figure; hold on; grid on; box on;
            h1 = plot( obs.t, Eus(inds,:), '-k' );
            h2 = plot( obs.t, EY(inds,:), '-r' );
            h3 = plot( obs.t, EP(inds,:), '-b' );
            h4 = plot( obs.t, Ediff(inds,:), ':r' );
            legend( [h1(1), h2(1), h3(1), h4(1)], ...
                {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
                 'Analytical (Paasschens, 1997)', 'diffusion approximation'}, ...
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
                             'polarization', 'P', ...
                             'lambda', 0.002 );
            
            material = MaterialClass.preset(3);
            
            observation = struct('r', 0:0.1:20, ... % size of bins in space
                                 'azimuth', [-pi pi], ... 
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
            obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
            
            Ep = squeeze(obs.energyDensity(:,:,:,1));
            Es = squeeze(obs.energyDensity(:,:,:,2));

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto( geometry, source, material, observation, false );
            %EtotY = sum(EY,3);

            % computing semi-analytical solution (Nakahara 2011, Sato 1994)
            % Eanalytical = Comparison.analyticalEnergyIsotropicElastic(geometry, material, observation, observation.r(inds), Wsp);
            
            % comparison of P energy densities
            figure; hold on; grid on; box on;
            h1 = plot( obs.t, Ep(inds,:), '-r' );
            h2 = plot( obs.t, EY(inds,:,1), '-b' );
            legend( [h1(1), h2(1)], ...
            {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'}, 'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('P-wave energy densities at different source-station distances')
            title(titlecase);

            % comparison of S energy densities
            figure; hold on; grid on; box on;
            h1 = plot( obs.t, Es(inds,:), '-r' );
            h2 = plot( obs.t, EY(inds,:,2), '-b' );
            legend( [h1(1), h2(1)], ...
            {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'}, 'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('S-wave energy densities at different source-station distances')
            title(titlecase);

        %% 3D Anisotropic scattering (isotropic differential scattering cross-section)
        % to be done
        case '3disotropicelastic'
            titlecase = '3D elastic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);

            geometry = struct( 'dimension', 3 );
            
            source = struct( 'numberParticles', 1e6, ...
                             'polarization', 'P', ...
                             'lambda', 0.002 );
            
            material = MaterialClass.preset(3);
            
            observation = struct('r', 0:0.1:20, ... % size of bins in space
                                 'azimuth', [-pi pi], ... 
                                 'elevation', [-pi/2 pi/2], ... 
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
            obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
            
            Ep = squeeze(obs.energyDensity(:,:,:,1))/obs.dphi;
            Es = squeeze(obs.energyDensity(:,:,:,2))/obs.dphi;
            Etotus = Ep + Es;

            % running Yoshimoto's Monte Carlo-based approach
            EY = Comparison.randomWalkYoshimoto( geometry, source, material, observation, false );
            EtotY = sum(EY,3);

            % computing semi-analytical solution (Nakahara 2011, Sato 1994)
            Eanalytical = Comparison.analyticalEnergyIsotropicElastic(geometry, material, observation, observation.r(inds), Wsp);
            
            % comparison of P energy densities
            figure; hold on; grid on; box on;
            h1 = plot( obs.t, Ep(inds,:), '-k' );
            h2 = plot( obs.t, EY(inds,:,1), '-b' );
            legend( [h1(1), h2(1)], ...
            {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'}, 'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('P-wave energy densities at different source-station distances')
            title(titlecase);

            % comparison of S energy densities
            figure; hold on; grid on; box on;
            h1 = plot( obs.t, Es(inds,:), '-k' );
            h2 = plot( obs.t, EY(inds,:,2), '-b' );
            legend( [h1(1), h2(1)], ...
            {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)'}, 'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('S-wave energy densities at different source-station distances')
            title(titlecase);

            % comparison of total energy densities
            figure; hold on; grid on; box on;
            h1 = plot( obs.t, Etotus(inds,:), '-k' );
            h2 = plot( obs.t, EtotY(inds,:), '-r' );
            h3 = plot( obs.t, Eanalytical, '-b' );
            legend( [h1(1), h2(1), h3(1)], ...
            {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
             'Analytical (Sato, 1994)'}, 'FontSize',12);
            xlabel('Lapse Time [s]');
            ylabel('Energy densities at different source-station distances')
            title(titlecase);

        otherwise
            error('unknown validation case')
    end
end