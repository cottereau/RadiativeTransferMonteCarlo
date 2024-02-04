function mainValidation(type)
% This function launches validation cases for RadiativeTransferMonteCarlo
%
% syntax mainValidation or mainValidation(validationCase)
%
% when validationCase='all' or empty, all validation cases are launched
%
% list of possible validationCases with reference to literature
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
    mainValidation('all')
else
    type = lower(type);
    switch type
        case 'all'
            mainValidation('2dIsotropicAcoustic')
            mainValidation('3dIsotropicAcoustic')
            mainValidation('2dIsotropicElastic')

            %% 2D Isotropic scattering acoustic (isotropic differential scattering cross-section)
        case '2disotropicacoustic'
            titlecase = '2D acoustic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            % input data
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
        case '3disotropicacoustic'
            titlecase = '3D acoustic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);
            % input data
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
            obs = radiativeTransferUnbounded( geometry.dimension, source, ...
                                                   material, observation );
            Eus = (obs.Ei+obs.Ec).*obs.dr';

            % computing Paasschens solution
            [EP,Ediff] = Paasschens_RTE_Unbounded( source, material, ...
                                                   observation, geometry );
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
        case '2disotropicelastic'
            titlecase = '2D elastic case with isotropic scattering';
            disp(['Testing ' titlecase ' ...']);

            % Load the results of the paper
            load('Nakahara_Yoshimoto_2011_Fig4_left.mat')


            % input data
            geometry = struct( 'dimension', 2 );
            source = struct( 'numberParticles', 1e5, ...
                'polarization', 'S', ...
                'lambda', 0.001 );
            observation = struct('dr', 0.4, ...        % size of bins in space
                'time', 0:0.43:34.13, ...  % observation times
                'Ndir', 10 );         % number of bins for directions
            material = struct( 'acoustics', false, ...
                'vp', 6, ...
                'vs', 3.46);

            % Values taken from (Nakahara & Yoshimoto, 2011, Section 3.1)
            K = material.vp/material.vs;
            gpp = 0.05; gps = 0.05; gp = gpp + gps;
            gsp = gps/K; gss = (material.vp*gp - material.vs*gsp)/material.vs;
            material.sigma = {@(th) 1/2/pi*ones(size(th))*gpp*material.vp, @(th) 1/2/pi*ones(size(th))*gps*material.vp; ...
                @(th) 1/2/pi*ones(size(th))*gsp*material.vs, @(th) 1/2/pi*ones(size(th))*gss*material.vs};

            % Run our code
            obsS = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
            source.polarizaton = 'P';
            obsP = radiativeTransferUnbounded( geometry.dimension, source, material, observation );

            % Epc = obs.Ec(:,:,1); Epi = obs.Ei(:,:,1); Ep = Epc+Epi; Esc = obs.Ec(:,:,2); Esi = obs.Ei(:,:,2); Es = Esc+Esi;
            EpcP = obsP.Ec(:,:,1); EpiP = obsP.Ei(:,:,1); EpP = EpcP+EpiP; EscP = obsP.Ec(:,:,2); EsiP = obsP.Ei(:,:,2); EsP = EscP+EsiP;
            EpcS = obsS.Ec(:,:,1); EpiS = obsS.Ei(:,:,1); EpS = EpcS+EpiS; EscS = obsS.Ec(:,:,2); EsiS = obsS.Ei(:,:,2); EsS = EscS+EsiS;
            Ep = (1*EpP+1.5*K^5*EpS)/(1+1.5*K^5); Es = (1*EsP+1.5*K^5*EsS)/(1+1.5*K^5);
            % Normalized P-wave energy in terms of normalized time
            figure; semilogy(obsP.t*gp*material.vp,2*pi*Ep(2,:)/(1/(1+1.5*K^5)*(gp*material.vp)^2),'-b','linewidth',2);
            xlabel('Normalized Time [-]'); ylabel('Normalized P-wave energy [-]');
            xlim([0 6]); ylim([0.001 100]); hold on; plot(energy_P(:,1),energy_P(:,2),'-r','linewidth',2);
            yticks([0.001 0.01 0.1 1 10 100]); yticklabels({'0.001','0.01','0.1','1','10','100'});

            % Normalized P-wave energy in terms of normalized time
            figure; semilogy(obsP.t*gp*material.vp,2*pi*Es(2,:)/(1/(1+1.5*K^5)*(gp*material.vp)^2),'-b','linewidth',2);
            xlabel('Normalized Time [-]'); ylabel('Normalized S-wave energy [-]');
            xlim([0 6]); ylim([0.001 100]); hold on; plot(energy_S(:,1),energy_S(:,2),'-r','linewidth',2);
            yticks([0.001 0.01 0.1 1 10 100]); yticklabels({'0.001','0.01','0.1','1','10','100'});

            %% 3D Anisotropic scattering (isotropic differential scattering cross-section)
            % to be done
        case '3dIsotropicElastic'

        otherwise
            error('unknown validation case')
    end
end

