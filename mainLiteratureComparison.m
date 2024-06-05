function mainLiteratureComparison(type)
% This function launches comparision cases for RadiativeTransferMonteCarlo
%
% syntax mainComparisonLiterature or mainComparisonLiterature(validationCase)
%
% when type='all' or empty, all validation cases are launched
%
% list of possible comparision cases with reference to literature
%  - '2dIsotropicAcoustic' (1,2)
%  - '3dIsotropicAcoustic' (1,2)
%  - '2dIsotropicElastic'  (3)
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

% with no argument, launch all validation cases
if nargin==0
    mainLiteratureComparison('all')
else
    type = lower(type);
    switch type
        case 'all'
            mainLiteratureComparison('2dIsotropicAcoustic')
            mainLiteratureComparison('3dIsotropicAcoustic')
            mainLiteratureComparison('2dIsotropicElastic')
            %% 2D Isotropic scattering acoustic (isotropic differential scattering cross-section)
        case '2disotropicacoustic'
            titlecase = '2D acoustic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            % input data
            source      = struct( 'numberParticles', 1e6, ...
                'position', [-10 0 0], ...
                'lambda', 0.001 );
            
            material = MaterialClass.preset(1);

            observation = struct( 'dr', 0.01, ...
                'time', 0:0.1:15, ...
                'Ndir', 100 );
            geometry    = struct( 'type', 'fullspace', ...
                'size', [4 3 3], ...
                'dimension', 2 );
            inds = [40 80 120]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
            Eus = (obs.Ei+obs.Ec).*obs.dr';

            % computing Paasschens solution
            [EP,Ediff] = Comparison.Paasschens_RTE_Unbounded( source, material, observation, geometry );

            % running Hoshiba's Monte Carlo-based approach
            EH = Comparison.Hoshiba_RTE_Unbounded_MonteCarlo( obs.t, obs.r(inds), ...
                source, material, geometry, 20 );
            EH = EH.*(pi*obs.dr(inds));

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
        case '3disotropicacoustic'
            titlecase = '3D acoustic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            % input data
            source      = struct( 'numberParticles', 1e6, ...
                'position', [-10 0 0], ...
                'lambda', 0.001 );
            material = MaterialClass.preset(2);
            observation = struct( 'dr', 0.05, ...
                'time', 0:0.1:30, ...
                'Ndir', 100 );
            geometry    = struct( 'type', 'fullspace', ...
                'size', [4 3 3], ...
                'dimension', 3 );
            inds = [40 80 120]; % index of the desired observation points

            % running our code, Monte Carlo-based
            obs = radiativeTransferUnbounded( geometry.dimension, source, ...
                material, observation );
            Eus = (obs.Ei+obs.Ec).*obs.dr';

            % computing Paasschens solution
            [EP,Ediff] = Comparison.Paasschens_RTE_Unbounded( source, material, ...
                observation, geometry );
            % % running Hoshiba's Monte Carlo-based approach
            % EH = Comparison.Hoshiba_RTE_Unbounded_MonteCarlo( obs.t, obs.r(inds), ...
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
        case '2disotropicelastic'
            titlecase = '2D elastic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);

            % Load the results of the paper
            % load('Nakahara_Yoshimoto_2011_Fig4_left.mat')

            % input data : values taken from (H. Nakahara and K. Yoshimoto 2011)
            geometry = struct( 'dimension', 2 );
            source = struct( 'numberParticles', 1e6, ...
                'polarization', 'S', ...
                'lambda', 0.001 );
            material = MaterialClass.preset(3);
            vp = material.vp; vs = material.vs;
            K = vp/vs;
            Sigmapp = 0.05*vp; Sigmaps = 0.05*vp; Sigmap = Sigmapp + Sigmaps;
            Sigmasp = Sigmaps/K^2; Sigmass = Sigmap-Sigmasp;
            %Sigmas = Sigmass + Sigmasp;
            Sigma = {Sigmapp, Sigmaps; Sigmasp, Sigmass};

            % if (Sigmap ~= Sigmas)
            %     keyboard
            %     error('This case cannot be dealt with in the framework of H. Sato 1994!')
            % else
            %     eta = Sigmap;
            % end
            eta = Sigmap;

            material.sigma = {@(th) 1/2/pi*ones(size(th))*Sigmapp, @(th) 1/2/pi*ones(size(th))*Sigmaps; ...
                @(th) 1/2/pi*ones(size(th))*Sigmasp, @(th) 1/2/pi*ones(size(th))*Sigmass};

            observation = struct('r', 0:0.1:20, ... % size of bins in space
                                 'time', 0:0.0167:10, ...       % observation times
                                 'Ndir', 10 );               % number of bins for directions

            % Run our code
            obsS = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

            source.polarization = 'P';
            obsP = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

            EpcP = obsP.Ec(:,:,1); EpiP = obsP.Ei(:,:,1); EpP = EpcP+EpiP;
            EscP = obsP.Ec(:,:,2); EsiP = obsP.Ei(:,:,2); EsP = EscP+EsiP;

            EpcS = obsS.Ec(:,:,1); EpiS = obsS.Ei(:,:,1); EpS = EpcS+EpiS;
            EscS = obsS.Ec(:,:,2); EsiS = obsS.Ei(:,:,2); EsS = EscS+EsiS;

            Ep = (1*EpP+1.5*K^5*EpS)/(1+1.5*K^5);
            Es = (1*EsP+1.5*K^5*EsS)/(1+1.5*K^5);

            % P-wave energy in terms of time
            b = 1; % Normalized distance
            ind = find(abs(obsP.r*eta/vp-b)<0.005);
            E_analytical = Comparison.analyticalEnergyIsotropicElastic(geometry.dimension,K,obsP.r(ind(1))*eta/vp,obsP.t,Sigma);
            figure; semilogy(eta*obsP.t,2*pi*Ep(ind(1),:)/(eta/vp)^2,'-b','linewidth',2);
            hold on; semilogy(eta*obsP.t,E_analytical(:,1),'-r','linewidth',2);
            xlabel('Time [s]'); ylabel('P-wave energy density');
            xlim([0 6]); ylim([1e-3 100]); grid on; box on;

            % S-wave energy in terms of time
            figure; semilogy(eta*obsP.t,2*pi*Es(ind(1),:)/(eta/vp)^2,'-b','linewidth',2);
            hold on; semilogy(eta*obsP.t,E_analytical(:,2),'-r','linewidth',2);
            xlabel('Time [s]'); ylabel('S-wave energy density');
            xlim([0 6]); ylim([1e-3 100]); grid on; box on;

            % Total energy in terms of time
            figure; plot(obsP.t,2*pi*(Ep(ind(1),:)+Es(ind(1),:))/(eta/vp)^2,'-b','linewidth',2);
            hold on; plot(obsP.t,sum(E_analytical,2),'-r','linewidth',2);
            xlabel('Time [s]'); ylabel('S-wave energy density');
            xlim([0 6]); grid on; box on;
            title(titlecase);

            %% 3D Anisotropic scattering (isotropic differential scattering cross-section)
            % to be done
        case '3disotropicelastic'
            titlecase = '3D elastic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);

            % input data : values taken from (H. Nakahara and K. Yoshimoto 2011)
            geometry = struct( 'dimension', 3 );
            source = struct( 'numberParticles', 1e6, ...
                'polarization', 'S', ...
                'lambda', 0.001 );
            material = MaterialClass.preset(4);
            vp = material.vp; vs = material.vs;
            K = vp/vs;
            Sigmapp = 0.05*vp; Sigmaps = 0.05*vp; Sigmap = Sigmapp + Sigmaps;
            Sigmasp = Sigmaps/K^2; Sigmass = Sigmap -Sigmasp;
            Sigmas = Sigmass + Sigmasp;
            Sigma = {Sigmapp, Sigmaps; Sigmasp, Sigmass};
            if (Sigmap ~= Sigmas)
                error('This case cannot be dealt with in the framework of H. Sato 1994!')
            else
                eta = Sigmap;
            end

            material.sigma = {@(th) 1/4/pi*ones(size(th))*Sigmapp, @(th) 1/4/pi*ones(size(th))*Sigmaps; ...
                @(th) 1/4/pi*ones(size(th))*Sigmasp, @(th) 1/4/pi*ones(size(th))*Sigmass};

            observation = struct('r', 0:0.05:15, ... % size of bins in space
                'time', 0:0.01:6, ...       % observation times
                'Ndir', 10 );               % number of bins for directions

            % Run our code
            obsS = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

            source.polarization = 'P';
            obsP = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

            EpcP = obsP.Ec(:,:,1); EpiP = obsP.Ei(:,:,1); EpP = EpcP+EpiP;
            EscP = obsP.Ec(:,:,2); EsiP = obsP.Ei(:,:,2); EsP = EscP+EsiP;

            EpcS = obsS.Ec(:,:,1); EpiS = obsS.Ei(:,:,1); EpS = EpcS+EpiS;
            EscS = obsS.Ec(:,:,2); EsiS = obsS.Ei(:,:,2); EsS = EscS+EsiS;

            Ep = (1*EpP+1.5*K^5*EpS)/(1+1.5*K^5);
            Es = (1*EsP+1.5*K^5*EsS)/(1+1.5*K^5);
            E = Ep + Es;

            % Total energy density in terms of time

            E_analytical = analyticalEnergyIsotropicElastic(geometry.dimension,K,[0.5 1 2],obsP.t,Sigma);
            figure; plot(obsP.t,2*pi*E(1,:),'-b','linewidth',2);
            hold on; plot(obsP.t,E_analytical(1,:)*(eta/vp)^3,'-r','linewidth',2);
            xlabel('Time [s]'); ylabel('Total energy density');
            xlim([0 6]); grid on; box on;
            title(titlecase);
        otherwise
            error('unknown validation case')
    end
end