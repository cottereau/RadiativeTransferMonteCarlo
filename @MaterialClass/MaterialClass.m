%% MaterialClass
% Store the physical, statistical, scattering, and transport properties of
% an acoustic or elastic medium. The class calculates the differential
% scattering cross-sections (DSCSs) required by the radiative transfer
% Monte Carlo solver for material parameter fluctuations and polycrystals.
classdef MaterialClass < handle
    properties
        % Wave physics and material description
        d                    int8   = 3;          % spatial dimension
        scatteringModel      char   = 'parameterFluctuations' % DSCS model
        correlationStructure char   = 'isotropic' % spatial correlation structure
        acoustics                                 % true: acoustic; false: elastic

        v                                         % acoustic-wave velocity [m/s]
        vp                                        % P-wave velocity [m/s]
        vs                                        % S-wave velocity [m/s]
        rho                                       % mass density [kg/m^3]
        Frequency                                 % carrier frequency [Hz]
        Q = Inf;                                  % quality factor(s); Inf means no intrinsic attenuation

        % Coefficients of variation: [kappa rho] for acoustics or
        % [lambda mu rho] for elasticity (dimensionless).
        coefficients_of_variation

        % Cross-correlation coefficients: corr(kappa,rho) for acoustics or
        % [corr(lambda,mu) corr(lambda,rho) corr(mu,rho)] for elastic waves.
        correlation_coefficients

        SpectralLaw          char   = '';          % PSDF/correlation-model name
        SpectralParam        struct = struct.empty % model-specific PSDF parameters
        CorrelationLength           = [];          % correlation length [m]

        % Polycrystal description
        singleCrystal        struct = struct.empty % single-crystal stiffness matrix and density
        TPCF                 struct = struct.empty % polycrystal two-point correlation function and spectrum

        % Scattering and transport properties
        sigma                cell   = cell.empty; % differential scattering cross-section(s) [1/s]
        Sigma                                     % total scattering cross-section(s) [1/s]
        Sigmapr                                   % forward-weighted cross-section(s) [1/s]
        invcdf                                    % inverse scattering-angle CDF(s)
        Diffusivity             double = [];      % diffusivity [m^2/s]
        meanFreeTime            double = [];      % mean free time(s) [s]
        meanFreePath            double = [];      % mean free path(s) [m]
        transportMeanFreeTime   double = [];      % transport mean free time(s) [s]
        transportMeanFreePath   double = [];      % transport mean free path(s) [m]
        g                       double = [];      % scattering anisotropy factor (dimensionless)
        P2P                                       % P-to-P scattering probability
        S2S                                       % S-to-S scattering probability

        % Spatial and spectral correlation functions
        Phi                           = [];       % normalized PSDF or function handle
        k                     double  = [];       % wavenumber sampling vector
        R                             = [];       % correlation function or function handle
        r                     double  = [];       % lag-distance sampling vector

        % Propagation-algorithm selector: 0 = small time steps; 1 = large time steps
        timeSteps                     = 0;

    end
    properties (Access = private)
        % Cached Zoeppritz amplitude and energy coefficients versus incidence angle
        zoeppritzOutputCache = [];
        % Cached functions giving reflected and transmitted wave angles
        zoeppritzAnglesCache = [];
        % Material properties used to determine whether the cache is still valid
        zoeppritzCacheKey    = [];
    end
    properties (SetAccess = private, Hidden = true)
        % Valid values accepted by the scatteringModel property set method
        ScatteringModel_def = {'parameterFluctuations','polycrystal'};
        % Valid values accepted by the correlationStructure property set method
        CorrelationStructure_def = {'isotropic','anisotropic'};
        % Valid values accepted by the SpectralLaw property set method
        SpectralLaw_def = {'','exp','power_law','gaussian','triangular', ...
            'low_pass','VonKarman','monodispersesphere','image','Imported'};
    end
    methods
        function obj = MaterialClass(geometry,freq,acoustics, ...
                v,coefficients_of_variation,correlation_coefficients, ...
                acf,lc)
            %% MaterialClass
            % Construct an acoustic or elastic material whose scattering is
            % caused by random fluctuations of its material parameters.
            %
            % Syntax:
            %   obj = MaterialClass()
            %   obj = MaterialClass(geometry,freq,acoustics,v, ...
            %       coefficients_of_variation,correlation_coefficients,acf,lc)
            %
            % Inputs:
            %   geometry : structure containing geometry.dimension
            %   freq     : carrier frequency [Hz]
            %   acoustics: true for acoustic waves; false for elastic waves
            %   v        : acoustic velocity or elastic velocities [Vp Vs] [m/s]
            %   coefficients_of_variation : [kappa rho] for acoustics or
            %                 [lambda mu rho] for elasticity
            %   correlation_coefficients : corr(kappa,rho) for acoustics or
            %                 [corr(lambda,mu) corr(lambda,rho) corr(mu,rho)]
            %   acf      : PSDF/correlation-model name
            %   lc       : correlation length [m]
            %
            % Output:
            %   obj      : configured MaterialClass object

            if nargin ~=0
                if geometry.dimension==1
                    warning(['Waves in 1D random media are localized and ' ...
                             'the RTT is not valid in this regime'])
                end

                obj.d = geometry.dimension;
                obj.acoustics = acoustics;
                obj.Frequency = freq;

                if acoustics
                    obj.v = v;
                else
                    obj.vp = v(1);
                    obj.vs = v(2);
                end

                obj.coefficients_of_variation = coefficients_of_variation;
                obj.correlation_coefficients = correlation_coefficients;

                obj.CorrelationLength = lc;
                obj.SpectralLaw = acf;

                obj.timeSteps = 0;
            end
        end
        function newobj = copyobj(obj)
            %% copyobj
            % Return a new class instance containing the same properties
            % values from the input
            %
            % Syntax
            %   newobj = copyobj(obj)
            %
            % Inputs:
            %   obj : Object to be copied
            %
            % Outputs:
            %   newobj : A new object from the same class with a copy of
            %   all its properties
            %
            if isscalar(obj)
                newobj = eval(mfilename);
                props = properties(obj);
                for i_props = 1:numel(props)
                    if isobject(obj.(props{i_props}))
                        if ~isempty(obj.(props{i_props}))
                            newobj.(props{i_props}) = obj.(props{i_props}).copyobj;
                        else
                            %empty
                        end
                    else
                        newobj.(props{i_props}) = obj.(props{i_props});
                    end
                end
            else
                for i=1:numel(obj)
                    newobj(i) = obj(i).copyobj; %#ok<AGROW>
                end
                newobj = reshape(newobj,size(obj));
            end
        end
        %% PROPERTY SET METHODS
        % MATLAB calls these methods automatically when the corresponding
        % property is assigned. Each method validates the requested value
        % and stores its canonical spelling.
        function set.scatteringModel(obj,newvalue)
            obj.scatteringModel = validatestring( ...
                newvalue,obj.ScatteringModel_def); %#ok<MCSUP>
        end
        function set.correlationStructure(obj,newvalue)
            obj.correlationStructure = validatestring( ...
                newvalue,obj.CorrelationStructure_def); %#ok<MCSUP>
        end
        function set.SpectralLaw(obj,newvalue)
            obj.SpectralLaw = validatestring(newvalue,obj.SpectralLaw_def); %#ok<MCSUP>
        end
        function CalcSigma(obj)
            %% CalcSigma
            % Compute the normalized differential scattering cross-sections based on
            % statistics of the fluctuating material parameters of the wave equation
            %
            % Syntax:
            %   newobj = CalcSigma (  );
            %
            % Inputs:
            %
            % Output:

            % The formulas are based on
            % L. Rhyzik, G. Papanicolaou, J. B. Keller. Transport equations for elastic
            % and other waves in random media. Wave Motion 24, pp. 327-370, 1996.
            % doi: 10.1016/S0165-2125(96)00021-2

            switch obj.scatteringModel
                case 'parameterFluctuations'
                    % Get the power spectral density function used by the
                    % generic material parameter fluctuation model.
                    if isempty(obj.Phi)
                        obj.getPSDF;
                    end

                    switch obj.correlationStructure
                        case 'isotropic'
                            obj.calcSigmaIsotropicCorrelation;
                        case 'anisotropic'
                            obj.calcSigmaAnisotropicCorrelation;
                        otherwise
                            error('Correlation structure not supported.')
                    end

                case 'polycrystal'
                    obj.calcSigmaPolycrystal;

                otherwise
                    error('Scattering model not supported.')
            end
        end
        function calcSigmaIsotropicCorrelation(obj)
            %% calcSigmaIsotropicCorrelation
            % Compute the normalized differential scattering cross-sections
            % based on the statistics of the fluctuating material parameters
            % of the wave equation for isotropic PSDF
            %
            % Syntax:
            %   newobj = calcSigmaIsotropicCorrelation( );
            %
            % Inputs:
            %
            % Output:
            % The 3D formulas are based on:
            %   L. Ryzhik, G. Papanicolaou, J. B. Keller. Transport equations for elastic
            %   and other waves in random media. Wave Motion 24, pp. 327-370, 1996.
            %   doi: 10.1016/S0165-2125(96)00021-2
            %
            % Note:
            % In Ryzhik et al, the power spectral density functions for 
            % lambda and mu are related to their respective reciprocals.
            % In our code, we consider lambda and mu fluctuations instead.
            % The cross PSDFs containing rho will thus change sign compared 
            % to the original formulas in the reference.
            %
            % The 2D elastic formulas correspond to the P-SV reduction of 
            % Ryzhik et al. elastic RTE and are consistent with the 
            % Born-scattering expressions in Sato et al. (2012), Chapter 4.

            switch obj.acoustics
                % Acoustic waves
                case 1
                    
                    omega = 2*pi*obj.Frequency;
                    zeta = omega/obj.v*obj.CorrelationLength;
                    delta_kk = obj.coefficients_of_variation(1); % CV of kappa (compressibility)
                    delta_rr = obj.coefficients_of_variation(2); % CV of rho (density)
                    rho_kr = obj.correlation_coefficients; % corr(kappa,rho)

                    switch obj.d
                        case 1
                            error('Waves in 1D random media are localized!')
                        case 2
                            % 2D ACOUSTICS
                            %
                            %   sigma(th) = (pi/2) * zeta^2 * omega
                            %               * (cos(th)^2 * delta_rr^2
                            %                  + 2*cos(th)*rho_kr*delta_kk*delta_rr
                            %                  + delta_kk^2)
                            %               * Phi(zeta * sqrt(2*(1-cos(th))))
                            %
                            % This is the 2D analogue of Ryzhik et al. (1996) Eq. (1.3).
                            % Note: Phi is the 2D isotropic PSDF evaluated at the
                            % scalar wavenumber |k' - k| = k*sqrt(2*(1-cos(th))).

                            obj.sigma = {@(th) (pi/2)*omega*zeta^2 * ...
                                (cos(th).^2*delta_rr^2 + ...
                                2*cos(th)*rho_kr*delta_kk*delta_rr + delta_kk^2) ...
                                .*obj.Phi(zeta.*sqrt(2*(1-cos(th))))};

                        case 3
                            % 3D ACOUSTICS
                            % Ryzhik et al, Eq. (1.3)
                            obj.sigma = {@(th) (pi/2)*omega*zeta^3*(cos(th).^2*delta_rr^2 + ...
                                2*cos(th)*rho_kr*delta_kk*delta_rr + delta_kk^2) ...
                                .*obj.Phi(zeta.*sqrt(2*(1-cos(th))))};
                    end
                
                % Elastic waves
                case 0

                    % Correlation inputs are direct physical fractional 
                    % fluctuations
                    K      = obj.vp / obj.vs;
                    omega  = 2*pi*obj.Frequency;
                    zetaP  = omega / obj.vp * obj.CorrelationLength;
                    zetaS  = K * zetaP;
                
                    delta_ll = obj.coefficients_of_variation(1); % CV of lambda
                    delta_mm = obj.coefficients_of_variation(2); % CV of mu
                    delta_rr = obj.coefficients_of_variation(3); % CV of rho
                
                    rho_lm   = obj.correlation_coefficients(1);  % corr(lambda,mu)
                    rho_lr   = obj.correlation_coefficients(2);  % corr(lambda,rho)
                    rho_mr   = obj.correlation_coefficients(3);  % corr(mu,rho)

                    switch obj.d
                        case 1
                            error('Waves in 1D random media are localized')
                        case 2
                        % 2D P-SV Elastics
                        %
                        % Note:
                        % Phi must be the 2D Ryzhik-normalized spectrum:
                        %   Rhat_2D(q) = variance * Lc^2 * Phi(q*Lc)
                        %
                        % If using Sato's raw 2D PSDF P_2D(q), where
                        %   P_2D(q) = variance * Lc^2 * Phi_Sato(q*Lc),
                        % then all prefactors below must be divided by (2*pi)^2.
                        %
                        % sigma{1,1}: P -> P
                        % sigma{1,2}: P -> SV
                        % sigma{2,1}: SV -> P
                        % sigma{2,2}: SV -> SV
                    
                        % ---------- P -> P ----------
                        sigmaPP = @(th) (pi/2)*omega*zetaP^2 .* ...
                            ( (1 - 2/K^2)^2*delta_ll^2 ...
                            + 4*(1/K^2 - 2/K^4)*rho_lm*delta_ll*delta_mm.*cos(th).^2 ...
                            + (4/K^4)*delta_mm^2.*cos(th).^4 ...
                            + delta_rr^2.*cos(th).^2 ...
                            - 2*(1 - 2/K^2)*rho_lr*delta_ll*delta_rr.*cos(th) ...
                            - (4/K^2)*rho_mr*delta_mm*delta_rr.*cos(th).^3 ) ...
                            .* obj.Phi(zetaP.*sqrt(2*(1-cos(th))));
                    
                        % ---------- P -> SV ----------
                        sigmaPS = @(th) (pi/2)*omega*zetaP^2 .* ...
                            ( K^2*delta_rr^2 ...
                            + 4*delta_mm^2.*cos(th).^2 ...
                            - 4*K*rho_mr*delta_mm*delta_rr.*cos(th) ) ...
                            .* (1-cos(th).^2) ...
                            .* obj.Phi(zetaP.*sqrt(1 + K^2 - 2*K*cos(th)));
                    
                        % ---------- SV -> P ----------
                        sigmaSP = @(th) (pi/2)*omega*zetaP^2 .* ...
                            ( delta_rr^2 ...
                            + (4/K^2)*delta_mm^2.*cos(th).^2 ...
                            - (4/K)*rho_mr*delta_mm*delta_rr.*cos(th) ) ...
                            .* (1-cos(th).^2) ...
                            .* obj.Phi(zetaS.*sqrt(1 + 1/K^2 - 2*cos(th)/K));
                    
                        % ---------- SV -> SV ----------
                        GammaSV = @(th) 2*cos(th).^2 - 1;
                    
                        sigmaSS = @(th) (pi/2)*omega*zetaS^2 .* ...
                            ( delta_rr^2.*cos(th).^2 ...
                            + delta_mm^2.*GammaSV(th).^2 ...
                            - 2*rho_mr*delta_mm*delta_rr.*cos(th).*GammaSV(th) ) ...
                            .* obj.Phi(zetaS.*sqrt(2*(1-cos(th))));
                    
                        obj.sigma = {sigmaPP, sigmaPS; ...
                                     sigmaSP, sigmaSS};
                        case 3
                            % 3D Elastics
                            %
                            % Note:
                            % Phi must be the 3D Ryzhik-normalized spectrum:
                            %   Rhat_3D(q) = variance * Lc^3 * Phi(q*Lc)
                            %
                            % If using Sato's raw 3D PSDF P_3D(q), where
                            %   P_3D(q) = variance * Lc^3 * Phi_Sato(q*Lc),
                            % then all prefactors below must be divided by (2*pi)^3.
                            %
                            % sigma{1,1}: P -> P
                            % sigma{1,2}: P -> S
                            % sigma{2,1}: S -> P
                            % sigma{2,2}: S -> S

                            % [Ryzhik et al, Eq. (1.3)] and [Turner, 1998, Eq. (3)]
                            sigmaPP = @(th) (pi/2)*omega*zetaP^3* ...
                                ( (1-2/K^2)^2*delta_ll^2 + 4*(1/K^2-2/K^4)*rho_lm*delta_ll*delta_mm*cos(th).^2 ...
                                + (4/K^4)*delta_mm^2*cos(th).^4 + delta_rr^2*cos(th).^2 ...
                                - 2*(1-2/K^2)*rho_lr*delta_ll*delta_rr*cos(th) ...
                                - (4/K^2)*rho_mr*delta_mm*delta_rr*cos(th).^3 ) ...
                                .*obj.Phi(zetaP.*sqrt(2*(1-cos(th))));

                            %  Ryzhik et al, Eqs. (4.56), (1.20), (1.22)
                            sigmaPS = @(th) (pi/2)*K*omega*zetaP^3* ...
                                ( K^2*delta_rr^2 + 4*delta_mm^2*cos(th).^2 - 4*K*rho_mr*delta_mm*delta_rr*cos(th) )...
                                .*(1-cos(th).^2).*obj.Phi(zetaP.*sqrt(1+K^2-2*K*cos(th)));

                            %  Ryzhik et al, Eqs. (4.56), (1.20), (1.21)
                            sigmaSP = @(th) (pi/4/K^3)*omega*zetaS^3* ...
                                ( delta_rr^2 + (4/K^2)*delta_mm^2*cos(th).^2 - (4/K)*rho_mr*delta_mm*delta_rr*cos(th) ) ...
                                .*(1-cos(th).^2).*obj.Phi(zetaS.*sqrt(1+1/K^2-2/K*cos(th)));

                            % Ryzhik et al, Eq. (4.54)
                            sigmaSS_TT = @(th) (pi/4)*omega*zetaS^3*delta_rr^2*(1+cos(th).^2)...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th))));
                            sigmaSS_GG = @(th) (pi/4)*omega*zetaS^3*delta_mm^2*(4*cos(th).^4-3*cos(th).^2+1)...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th))));
                            sigmaSS_GT = @(th) -(pi/4)*omega*zetaS^3*rho_mr*delta_mm*delta_rr*(4*cos(th).^3)...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th))));

                            sigmaSS = @(th) sigmaSS_TT(th) + sigmaSS_GG(th) + sigmaSS_GT(th);

                            obj.sigma = {sigmaPP,sigmaPS; ...
                                         sigmaSP,sigmaSS};
                    end
            end
        end
        function calcSigmaAnisotropicCorrelation(~)
            %% calcSigmaAnisotropicCorrelation
            % Compute the normalized differential scattering cross-sections
            % based on statistics of the fluctuating material parameters
            % of the wave equation for anisotropic PSDF
            %
            % Syntax:
            %   newobj = calcSigmaAnisotropicCorrelation (  );
            %
            % Inputs:
            %
            % Output:

            % The formulas are based on
            % L. Rhyzik, G. Papanicolaou, J. B. Keller. Transport equations for elastic
            % and other waves in random media. Wave Motion 24, pp. 327-370, 1996.
            % doi: 10.1016/S0165-2125(96)00021-2


            error('Not implemented yet')
        end
        function calcSigmaPolycrystal(obj)
            %% calcSigmaPolycrystal
            % Calculate the four 3D polycrystal differential scattering
            % cross-sections required by the radiative-transfer solver.
            %
            % The spatial cross-sections (Omega) follow the generalized 
            % Weaver framework. Multiplication by the incident phase velocity 
            % gives the time-based operators stored in obj.sigma. 
            % MaterialClass.prepareSigma subsequently calculates
            % Sigma, inverse angular CDFs, and the mean free quantities.

            if obj.d ~= 3
                error('MaterialClass:PolycrystalDimension', ...
                      ['The polycrystal scattering model currently requires ', ...
                       'd = 3.']);
            end
            if isempty(obj.acoustics) || obj.acoustics
                error('MaterialClass:PolycrystalElasticOnly', ...
                      'The polycrystal scattering model requires elastic waves.');
            end
            validateattributes(obj.Frequency, {'numeric'}, ...
                {'real','finite','scalar','positive'}, ...
                 'MaterialClass.calcSigmaPolycrystal', 'Frequency');

            if ~isscalar(obj.singleCrystal) || ...
               ~isfield(obj.singleCrystal,'rho') || ...
               ~isfield(obj.singleCrystal,'C')
                  error('MaterialClass:InvalidSingleCrystal', ...
                        'singleCrystal must be a scalar structure containing rho and C.');
            end

            density = obj.singleCrystal.rho;
            Cloc = obj.singleCrystal.C;
            validateattributes(density, {'numeric'}, ...
                {'real','finite','scalar','positive'}, ...
                'MaterialClass.calcSigmaPolycrystal', 'singleCrystal.rho');
            if ~isnumeric(Cloc) || ~isreal(Cloc) || ...
                    ~isequal(size(Cloc),[6 6]) || any(~isfinite(Cloc),'all')
                error('MaterialClass:InvalidSingleCrystalStiffness', ...
                    'singleCrystal.C must be a finite, real 6-by-6 matrix.');
            end
            if ~isscalar(obj.TPCF) || ~isfield(obj.TPCF,'spectrum') || ...
                    ~isa(obj.TPCF.spectrum,'function_handle')
                error('MaterialClass:InvalidPolycrystalTPCF', ...
                    'TPCF.spectrum must be a function handle of wavenumber q.');
            end

            % Corrected Kube-Turner directional inner products
            [L, M, N] = MaterialClass.InnerProducts(Cloc);

            % Calculate Voigt velocities unless the user prescribed them
            if isempty(obj.vp) && isempty(obj.vs)
                C = 0.5*(Cloc + Cloc.');
                normalSum = C(1,1) + C(2,2) + C(3,3);
                normalCrossSum = C(1,2) + C(1,3) + C(2,3);
                shearSum = C(4,4) + C(5,5) + C(6,6);
                cP = (3*normalSum + 2*normalCrossSum + 4*shearSum)/15;
                cS = (normalSum - normalCrossSum + 3*shearSum)/15;
                if cP <= 0 || cS <= 0
                    error('MaterialClass:NonPositiveVoigtModulus', ...
                          'The Voigt P- and S-wave moduli must be positive.');
                end
                obj.vp = sqrt(cP/density);
                obj.vs = sqrt(cS/density);
            else
                velocities = [obj.vp obj.vs];
                if numel(velocities) ~= 2 || any(~isfinite(velocities)) || ...
                   any(velocities <= 0)
                     error('MaterialClass:InvalidPolycrystalVelocities', ...
                           'Prescribed phase velocities must be positive [Vp Vs].');
                end
            end
            obj.rho = density;

            vP = obj.vp;
            vS = obj.vs;
            omega = 2*pi*obj.Frequency;
            kp = omega/vP;
            ks = omega/vS;
            etaTilde = obj.TPCF.spectrum;

            % GSD (grain size distribution) averages are more expensive 
            % than closed-form spectra.
            % Evaluate them once over the complete wavenumber-transfer
            % interval needed at this frequency, then use a linear lookup
            % while constructing and integrating the scattering operators.
            if isfield(obj.TPCF,'model') && ...
                    any(strcmp(obj.TPCF.model, ...
                    {'Arguelles','Sha','ShengKhazaie'}))
                % For q = |k_scattered - k_incident|, the maximum occurs at
                % backscattering (theta = pi), where q = k_incident+k_scattered.
                % Across PP, PS, SP, and SS scattering, qmax = 2*max(kp,ks).
                maximumTransferWavenumber = 2*max(kp,ks);
                transferWavenumber = linspace(0,maximumTransferWavenumber,4097);
                spectralValue = etaTilde(transferWavenumber);
                spectralInterpolant = griddedInterpolant(transferWavenumber,spectralValue,'linear','nearest');
                etaTilde = @(q) spectralInterpolant(abs(q));
            end

            mu = @(theta) cos(theta);
            IPpp = @(theta) L.L0 + L.L1.*mu(theta).^2 + L.L2.*mu(theta).^4;
            IPps = @(theta) M.M0 + M.M1.*mu(theta).^2 - IPpp(theta);
            IPsp = IPps;
            IPss = @(theta) N.N0 + N.N1.*mu(theta).^2 - 2*(M.M0 + M.M1.*mu(theta).^2) + IPpp(theta);

            qpp = @(theta) kp.*sqrt(max(0,2*(1-mu(theta))));
            qps = @(theta) sqrt(max(0,kp^2 + ks^2 - 2*kp*ks.*mu(theta)));
            qsp = qps;
            qss = @(theta) ks.*sqrt(max(0,2*(1-mu(theta))));

            prefactorPP = pi*omega^4/(4*density^2*vP^8);
            prefactorPS = pi*omega^4/(4*density^2*vP^3*vS^5);
            prefactorSP = pi*omega^4/(8*density^2*vS^3*vP^5);
            prefactorSS = pi*omega^4/(8*density^2*vS^8);

            OmegaPP = @(theta) prefactorPP.*IPpp(theta).*etaTilde(qpp(theta));
            OmegaPS = @(theta) prefactorPS.*IPps(theta).*etaTilde(qps(theta));
            OmegaSP = @(theta) prefactorSP.*IPsp(theta).*etaTilde(qsp(theta));
            OmegaSS = @(theta) prefactorSS.*IPss(theta).*etaTilde(qss(theta));

            obj.sigma = {@(theta) vP.*OmegaPP(theta), ...
                         @(theta) vP.*OmegaPS(theta); ...
                         @(theta) vS.*OmegaSP(theta), ...
                         @(theta) vS.*OmegaSS(theta)};

            % Check vectorization, finiteness, and physical non-negativity
            theta = linspace(0,pi,1001);
            for incident = 1:2
                for scattered = 1:2
                    value = obj.sigma{incident,scattered}(theta);
                    if ~isnumeric(value) || ~isreal(value) || ...
                            ~isequal(size(value),size(theta)) || ...
                            any(~isfinite(value))
                        error('MaterialClass:InvalidPolycrystalDSCS', ...
                             ['Each polycrystal DSCS must return a finite, ', ...
                              'real array with the same size as theta.']);
                    end
                    tolerance = 1e-10*max(1,max(abs(value)));
                    if any(value < -tolerance)
                        error('MaterialClass:NegativePolycrystalDSCS', ...
                            'A polycrystal differential cross-section is negative.');
                    end
                end
            end
        end
        function getPSDF(obj)
            % Compute the normalized power spectral density function to use
            % it in the differential scattering cross-section
            %
            % Syntax:
            %   getPSDF (  );
            %
            % Inputs:
            %
            % Output:

            % In Ryzhik et al, the PSDFs are those of the fractional parts 
            % (normalized) of the corresponding random fields.
            % We assume a factorizable spectral covariance:
            %   Rhat_ab(q) = sigma_a * sigma_b * rho_ab * Phi(q)
            % where sigma_a are the coefficients of variation and rho_ab are
            % correlation coefficients. Here we check only the correlation matrix
            % with ones on the diagonal and rho_ab off diagonal.
            % Extending to the non-factorizable case is straightforward.

            if obj.acoustics
                rho_kr = obj.correlation_coefficients;
            
                if abs(rho_kr) > 1
                    error('Absolute value of corr(kappa,rho) should be less than 1!')
                end
            
            else
                rho_lm = obj.correlation_coefficients(1);
                rho_lr = obj.correlation_coefficients(2);
                rho_mr = obj.correlation_coefficients(3);
            
                Ccorr = [1      rho_lm  rho_lr; ...
                         rho_lm 1       rho_mr; ...
                         rho_lr rho_mr  1     ];
            
                if min(eig(Ccorr)) < -1e-12
                    error('The lambda-mu-rho correlation matrix is not positive semidefinite.');
                end
            
                if any(abs(obj.correlation_coefficients)>1)
                    error('Absolute values of correlation coefficients should be less than 1!')
                end
            end

            % warning: these are only for 3D (formulas depend on the dimension)
                    switch obj.SpectralLaw
                        case 'exp'
                            obj.Exponential(obj.CorrelationLength);
                        case 'power_law'
                            obj.PowerLaw(obj.CorrelationLength);
                        case 'gaussian'
                            obj.Gaussian(obj.CorrelationLength);
                        case 'triangular'
                            obj.Triangular(obj.CorrelationLength);
                        case 'low_pass'
                            obj.LowPass(obj.CorrelationLength);
                        case 'VonKarman'
                            if ~isfield(obj.SpectralParam,'nu')
                                error('Please add the Hurst number for VonKarman PSDF')
                            end
                            obj.VonKarman(obj.CorrelationLength, obj.SpectralParam.nu);
                        case 'monodispersesphere'
                            if ~isfield(obj.SpectralParam,'rhoS') && ~isfield(obj.SpectralParam,'Diam')
                                error('Please add the mono disperse disk parameters for the PSDF')
                            end
                            obj.MonoDisperseSphere(obj.SpectralParam.rhoS,obj.SpectralParam.Diam);
                        case 'image'
                            if ~isfield(obj.SpectralParam,'ImagePath') && ~isfield(obj.SpectralParam,'dx') &&~isfield(obj.SpectralParam,'dy')
                                error('not defined yet')
                            end
                            obj.GetPSDFromImage(obj.SpectralParam);
                    end

        end

        % Correlation-length convention:
        %   1D: Lc   = 2 * int_0^inf R(r) dr
        %   2D: Lc^2 = 2 * int_0^inf r * R(r) dr
        %   3D: Lc^3 = 3 * int_0^inf r^2 * R(r) dr
        %
        % In the analytical models below, the dimensionless correlation
        % function R(z), z = r/Lc, should satisfy Lc = 1 under the
        % corresponding formula above.
        %
        % Fourier-transform convention:
        %   Phi(k) = 1/(2*pi)^d * int_{R^d} exp(-i k.x) R(|x|) dx.
        %
        % For isotropic correlation functions, this becomes:
        %   1D: Phi(k) = 2/(2*pi) * int_0^inf cos(k*r) R(r) dr
        %   2D: Phi(k) = 2*pi/(2*pi)^2 * int_0^inf r*J0(k*r) R(r) dr
        %   3D: Phi(k) = 4*pi/(2*pi)^3 * int_0^inf r^2*sinc(k*r) R(r) dr
        %
        % with sinc(x) = sin(x)/x.
        %
        % With this convention, the dimensional spectrum is
        %   Rhat_d(q) = sigma^2 * Lc^d * Phi(q*Lc).

        function out = Exponential(obj,lc)
            %% Exponential
            % Compute the normalized power spectral density function for
            % Exponential
            %
            % Syntax:
            %   Exponential (  );
            %
            % Inputs:
            %  lc: correlation length
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            obj.SpectralParam = struct('correlationLength',lc);
            obj.CorrelationLength = lc;
            obj.SpectralLaw = 'Exp';

            if obj.d == 1
                obj.R = @(z) exp(-2*z);
                obj.Phi = @(z) 1./(2*pi*(1+(z/2).^2));
            elseif obj.d == 2
                a = sqrt(2);
                obj.R = @(z) exp(-a*z);
                obj.Phi = @(z) 1./(4*pi*(1+(z/a).^2).^(1.5));
            
            elseif obj.d == 3
                a = 6^(1/3);
                obj.R = @(z) exp(-a*z);
                obj.Phi = @(z) 1./(6*pi^2*(1+(z/a).^2).^2);
            else
                error('incorrect dimension!')
            end
            %out = @(z) 1./(8*pi^2*(1+z.^2/4).^2);
            %obj.R = @(z) exp(-2*z);
            out = obj.Phi;
        end
        
        function out = PowerLaw(obj,lc)
            %% PowerLaw
            % Compute the normalized power spectral density function for
            % Power Law
            %
            % Syntax:
            %   PowerLaw (  );
            %
            % Inputs:
            %  lc: correlation length
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            obj.SpectralParam = struct('correlationLength',lc);
            obj.CorrelationLength = lc;
            obj.SpectralLaw = 'power_law';
            % out = @(z) 1./(pi^4)*exp(-2*z/pi);
            % obj.Phi = out;
            % obj.R = @(z) 1./(1+(pi^2*z.^2/4))^2;
            if obj.d == 1
                a = pi/2;
                obj.R = @(z) 1./(1+(a*z).^2).^2;
                obj.Phi = @(z) 1./(2*pi).*(1+z/a).*exp(-z/a);
            
            elseif obj.d == 2
                obj.R = @(z) 1./(1+z.^2).^2;
                obj.Phi = @(z) z.*besselk(1,z)./(4*pi);
                % use the limit Phi(0)=1/(4*pi) if needed
            
            elseif obj.d == 3
                a = (3*pi/4)^(1/3);
                obj.R = @(z) 1./(1+(a*z).^2).^2;
                obj.Phi = @(z) 1./(6*pi^2).*exp(-z/a);
            end

            out = obj.Phi;
        end

        function out = Gaussian(obj,lc)
            %% Gaussian
            % Compute the normalized power spectral density function for
            % Gaussian
            %
            % Syntax:
            %   Gaussian (  );
            %
            % Inputs:
            %  lc: correlation length
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            obj.SpectralParam = struct('correlationLength',lc);
            obj.CorrelationLength = lc;
            obj.SpectralLaw = 'Gauss';

            if obj.d == 1
                obj.R = @(z) exp(-pi*z.^2);
                obj.Phi = @(z) 1./(2*pi).*exp(-z.^2/(4*pi));
            elseif obj.d == 2
                obj.R = @(z) exp(-z.^2);
                obj.Phi = @(z) 1./(4*pi).*exp(-z.^2/4);
            elseif obj.d == 3
                a = (3*sqrt(pi)/4)^(2/3);
                obj.R = @(z) exp(-a*z.^2);
                obj.Phi = @(z) 1./(6*pi^2).*exp(-z.^2/(4*a));
            else
                error('incorrect dimension!')
            end

            out = obj.Phi;

            % out = @(z) 1./(8*pi^3)*exp(-z.^2/4/pi);
            % obj.Phi = out;
            % obj.R = @(z)exp(-pi*z.^2);
        end

        function out = Triangular(obj,lc)
            %% Triangular
            % Compute the normalized power spectral density function for
            % Triangular
            %
            % Syntax:
            %   Triangular (  );
            %
            % Inputs:
            %  lc: correlation length
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            obj.SpectralParam = struct('correlationLength',lc);
            obj.CorrelationLength = lc;
            obj.SpectralLaw = 'Triangular';

            if obj.d == 3
                obj.R = @(z)(12*(2-2*cos(2*pi*z)-(2*pi*z).*sin(2*pi*z)))./(2*pi*z).^4;
                obj.Phi = @(z) (3/8/pi^4)*(1-z/2/pi).*obj.heaviside(2*pi-z);
            else
                error('dimensions other than 3 are not coded yet!');
            end

            out = obj.Phi;
            
        end

        function out = LowPass(obj,lc)
            %% LowPass
            % Compute the normalized power spectral density function for
            % LowPass
            %
            % Syntax:
            %   LowPass (  );
            %
            % Inputs:
            %  lc: correlation length
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            obj.SpectralParam = struct('correlationLength',lc);
            obj.CorrelationLength = lc;
            obj.SpectralLaw = 'low_pass';
            
            if obj.d == 3
                a = 3*pi/2;
                obj.R = @(z) (3*(sin(a*z)-a*z.*cos(a*z)))./(a*z).^3;
                obj.Phi = @(z) (2/9/pi^4)*obj.heaviside(a-z);
            else
                error('LowPass spectral model is currently implemented only for d = 3.');
            end

            out = obj.Phi;

        end

        function out = VonKarman(obj,lc,nu)
            %% VonKarman
            % Compute the normalized power spectral density function for
            % Von Karman
            %
            % Syntax:
            %   VonKarman (  );
            %
            % Inputs:
            %   nu : Hurst number
            %
            % Output:

            % https://reproducibility.org/RSF/book/sep/fractal/paper_html/node4.html
            % Goff, J. A., and T. H. Jordan, 1988, Stochastic modeling of
            % seafloor morphology: Inversion of sea beam data for
            % second-order statistics: Journal of Geophysical Research,
            % 93, 13,589-13,608.

            % $K_\nu$ is the modified Bessel function of order $\nu $,
            % where $0.0<\nu<1.0$ is the Hurst number (Mandelbrot, 1985,1983).

            obj.SpectralParam = struct('correlationLength',lc,'nu',nu);
            obj.SpectralLaw = 'VonKarman';
            obj.CorrelationLength = lc;

            if obj.d == 1
                a = 2*sqrt(pi)*gamma(nu+0.5)/gamma(nu);
                obj.R = @(z) 2^(1-nu)/gamma(nu) * (a*z).^nu .* besselk(nu,a*z);
                obj.Phi = @(z) 1/(2*pi) ./ (1+(z/a).^2).^(nu+0.5);
            
            elseif obj.d == 2
                a = 2*sqrt(nu);
                obj.R = @(z) 2^(1-nu)/gamma(nu) * (a*z).^nu .* besselk(nu,a*z);
                obj.Phi = @(z) 1/(4*pi) ./ (1+(z/a).^2).^(nu+1);
            
            elseif obj.d == 3
                a = (6*sqrt(pi)*gamma(nu+1.5)/gamma(nu))^(1/3);
                obj.R = @(z) 2^(1-nu)/gamma(nu) * (a*z).^nu .* besselk(nu,a*z);
                obj.Phi = @(z) 1/(6*pi^2) ./ (1+(z/a).^2).^(nu+1.5);
            else
                error('incorrect dimension!')
            end

            out = obj.Phi;

        end

        function out = MonoDisperseSphere(obj, eta, D)
            %% MonoDisperseSphere
            % Compute the normalised power spectral density function for a
            % monodisperse hard-body suspension using Percus-Yevick theory.
            %
            % The pair correlation function g(r) and static structure factor S(k)
            % are computed by MaterialClass.HardBodyPY for d = 1, 2, or 3.
            % The two-point probability S2(r) and PSDF Phi(k) are then derived
            % from g(r)/S(k) following the Torquato spectral route.
            %
            % References:
            %   Torquato & Stell (1985) J. Chem. Phys. 82(2), 980-987  [3D]
            %   Adda-Bedia, Katzav & Vella (2008) J. Chem. Phys. 128, 184508  [2D]
            %   Torquato, Random Heterogeneous Materials, Springer (2002)
            %
            % Syntax:
            %   out = MonoDisperseSphere(obj, eta, D)
            %
            % Inputs:
            %   eta : packing fraction (length/area/volume fraction for d=1/2/3)
            %   D   : object diameter (same units as desired length scale)
            %
            % Output:
            %   out : function handle for the PSDF  Phi(k)  (normalised)

            switch obj.d
                case 1
                    error('Not implemented')

                case 2
                    % =========================================================
                    % 1. Parameters
                    % =========================================================
                    %D        = 1.0;          % disk diameter (unit length)
                    %phi      = 0.45;         % area packing fraction
                    % phi = eta;
                    rho      = eta / (pi*(D/2)^2);  % number density [disks / area]

                    Nr       = 4096;         % number of real-space grid points
                    r_max    = 20.0 * D;     % real-space cutoff
                    dr       = r_max / Nr;
                    r        = ((1:Nr) - 0.5)' * dr;   % cell-centred grid, avoids r=0

                    max_iter = 5000;         % max Picard iterations
                    tol      = 1e-12;        % convergence tolerance on max|delta_c|
                    alpha    = 0.4;          % Picard mixing parameter (0 < alpha <= 1)

                    fprintf('=========================================\n');
                    fprintf(' 2D Percus-Yevick OZ solver — hard disks\n');
                    fprintf('=========================================\n');
                    fprintf(' Diameter D    = %.4f\n', D);
                    fprintf(' Packing phi   = %.4f\n', eta);
                    fprintf(' Number density= %.4f\n', rho);
                    fprintf(' Grid points   = %d,  dr = %.5f\n', Nr, dr);

                    % =========================================================
                    % 2. Hard-disk pair potential  (beta*u = 0 outside, inf inside)
                    %    PY hard-core condition: g(r) = 0  for r < D
                    %    => c(r) = -1              for r < D  (exact PY hard-core)
                    %    => c(r) = h(r) - ln g(r)  for r > D  (PY closure, hard pot.)
                    %    Simplifies to: c(r) = g(r) - 1 - h(r)*[g(r)-1]/g(r) ...
                    %    Standard form used here:
                    %       c(r) = (1 + h(r)) * (1 - exp(beta*u)) = 0   for r > D (hard)
                    %    => c(r) = g(r) - 1 - [g(r)-1]     ...
                    %    Cleanest hard-disk PY:
                    %       for r <  D : g=0, c = -1  (from h = -1, PY => c = g*y - y = -y; y=exp(-bu)*g; for hard core outside: c=g-1-h*(g-1)/g ... )
                    %    We use the standard iterative form directly:
                    %       PY closure outside core: c(r) = g(r) - 1 - gamma(r)
                    %       where gamma(r) = h(r) - c(r)  is the indirect correlation fn
                    %       Inside core: g(r)=0 => h(r)=-1 => c(r) = -1 - gamma(r)
                    % =========================================================
                    inside  = r < D;   % logical mask: inside hard core
                    outside = ~inside;

                    % =========================================================
                    % 3. Build Hankel transform matrix / use discrete sine-like approach
                    %    For 2D isotropic functions, the Fourier transform is the
                    %    zeroth-order Hankel transform:
                    %       f_hat(k) = 2*pi * int_0^inf f(r) J0(k*r) r dr
                    %    On a uniform grid we use the quasi-discrete Hankel transform (QDHT)
                    %    or simply a direct quadrature via the trapezoidal rule with J0.
                    %    For efficiency we build a uniform k-grid matching the r-grid and
                    %    use the FFT-based approach via the identity:
                    %       H0{f}(k) = 2*pi * int f(r) J0(kr) r dr
                    %    Implemented here as a direct matrix-free quadrature using
                    %    the fast oscillatory integral approximation via FFT + Abel:
                    %       We convert to h_tilde(k) = 2*pi * DHTQ(r*f(r), k)
                    %    For simplicity and robustness we use direct Gauss quadrature
                    %    with the J0 kernel (accurate for smooth functions on [0, r_max]).
                    % =========================================================
                    %  k-grid: same spacing as r-grid for reciprocal consistency
                    dk    = pi / r_max;        % Nyquist-consistent spacing
                    k_max = pi / dr;
                    Nk    = floor(k_max / dk);
                    k     = ((1:Nk) - 0.5)' * dk;

                    fprintf(' k-grid: Nk=%d, dk=%.5f, k_max=%.2f\n', Nk, dk, k_max);

                    % Precompute Hankel kernel: H(i,j) = 2*pi * J0(k(i)*r(j)) * r(j) * dr
                    % This is memory-intensive for large Nr*Nk; we use chunked evaluation.
                    % For Nr=Nk=4096 the full matrix is 4096^2 * 8 bytes = 128 MB — manageable.
                    fprintf(' Building Hankel kernel ... ');
                    % We do NOT build the full matrix. Instead we define inline functions
                    % that apply the transform via a loop over k-chunks (memory efficient).

                    hankel_fwd = @(f_r) MaterialClass.hankel0_fwd(f_r, r, k, dr);
                    hankel_inv = @(f_k) MaterialClass.hankel0_fwd(f_k, k, r, dk) / (2*pi)^2 * (2*pi);
                    % Note: inverse Hankel in 2D = (1/(2*pi)) * H0, so:
                    %   f(r) = (1/(2*pi)) * int f_hat(k) J0(kr) k dk
                    hankel_inv = @(f_k) MaterialClass.hankel0_inv(f_k, k, r, dk);
                    fprintf('done\n');

                    % =========================================================
                    % 4. Initial guess: c(r) from low-density limit
                    %    c(r) ~ -1  for r < D, c(r) ~ 0 for r > D
                    % =========================================================
                    gamma_r = zeros(Nr, 1);   % indirect correlation: gamma = h - c
                    c_r     = zeros(Nr, 1);
                    c_r(inside)  = -1.0;      % hard core: exact PY value at r < D
                    c_r(outside) =  0.0;      % dilute start

                    fprintf('\n Starting Picard iteration (alpha=%.2f, tol=%.1e)...\n', alpha, tol);

                    % =========================================================
                    % 5. Picard iteration
                    % =========================================================
                    converged = false;
                    for iter = 1:max_iter

                        % --- OZ in k-space ---
                        c_k   = hankel_fwd(c_r);
                        h_k   = c_k ./ (1.0 - rho * c_k);   % OZ solution
                        h_r   = hankel_inv(h_k);

                        % --- Update gamma ---
                        gamma_r_new = h_r - c_r;

                        % --- PY closure for new c ---
                        c_r_new = zeros(Nr, 1);
                        % Inside core: g=0 => h=-1 => c = -1 - gamma
                        c_r_new(inside) = -1.0 - gamma_r(inside);
                        % Outside core: hard potential => exp(-beta*u)=1
                        %   PY: c = (g)(1 - exp(beta*u)) = 0  ... cleaner:
                        %   c = h - gamma*(g)  => for hard outside: c = (1+gamma)*(1) - 1 - gamma ...
                        %   Standard: c_out = (1+gamma)*f  where f=exp(-bu)-1=0 for hard outside => c_out = 0
                        %   But correct PY for hard disk outside: c(r) = g(r) - 1 - gamma(r)*...
                        %   Using h = gamma + c and PY => c = (h+1)*[exp(-bu)-1]/(...
                        %   Correct hard-disk PY outside core:
                        %     c(r) = -gamma(r)  * [g(r)=1+gamma(r)+c(r)]/(...)
                        %   Simplest consistent form (Lado 1968):
                        %     c(r) = (1 + gamma_r) - 1 - gamma_r = 0   ... NO
                        %   The exact PY closure is: c(r) = exp(-beta*u(r)) * (1 + gamma(r)) - 1 - gamma(r)
                        %   For hard disks outside (u=0): c(r) = (1+gamma) - 1 - gamma = 0 ???
                        %   That gives c=0 outside always which is wrong at higher density.
                        %   The issue: PY closure is  c(r) = g(r)*(1 - exp(+beta*u(r)))
                        %   For hard outside (u=0): c(r) = g(r)*0 = 0  -- this IS the PY result.
                        %   The indirect correlation gamma carries all the info outside the core.
                        %   So: outside core, c_new = 0 always in PY for hard disks.
                        c_r_new(outside) = 0.0;

                        % --- Convergence check ---
                        err = max(abs(gamma_r_new - gamma_r));
                        if mod(iter, 200) == 0
                            fprintf('  iter %4d : err = %.3e\n', iter, err);
                        end

                        % --- Picard mixing ---
                        gamma_r = (1 - alpha) * gamma_r + alpha * gamma_r_new;

                        % Reconstruct c from mixed gamma
                        c_r(inside)  = -1.0 - gamma_r(inside);
                        c_r(outside) = 0.0;

                        if err < tol
                            fprintf('  Converged at iter %d, err = %.3e\n', iter, err);
                            converged = true;
                            break;
                        end
                    end

                    if ~converged
                        fprintf('  WARNING: did not fully converge (err=%.3e). Using current solution.\n', err);
                    end

                    % =========================================================
                    % 6. Final correlation functions
                    % =========================================================
                    c_k   = hankel_fwd(c_r);
                    h_k   = c_k ./ (1.0 - rho * c_k);
                    h_r   = hankel_inv(h_k);
                    g_r   = 1.0 + h_r;            % pair correlation function (RDF)
                    g_r(inside) = 0.0;             % enforce hard core

                    % Structure factor S(k)
                    S_k   = 1.0 + rho * h_k;      % S(k) = 1 + rho * h_hat(k)

                    % =========================================================
                    % 7. Two-point correlation function S2(r)
                    %    For a stationary, isotropic 2D medium:
                    %    S2(r) = phi^2 + phi*(1-phi)^2 * [g(r)-1]   [dilute approx]
                    %    Exact form (Torquato 2002, Eq. 2.44):
                    %    S2(r) = phi^2 [1 + (1/phi^2) * rho^2 * int int ...]
                    %    In practice, for a single-phase indicator:
                    %    S2(r) = phi^2 + rho^2 * int h(r12) m(r1) m(r2) dr1 dr2 / V
                    %    where m is the single-disk indicator.
                    %    A common tractable approximation (valid at all densities):
                    %
                    %    chi(r) = S2(r) - phi^2  = phi*(1-phi)*[rho*h_hat convolved with m_hat^2 / V]
                    %
                    %    The "spectral density" approach (Torquato & Lu):
                    %    chi_hat(k) = rho * |m_hat(k)|^2 * S(k)
                    %    where m_hat(k) is the 2D Fourier transform of a single disk indicator:
                    %    m_hat(k) = pi*D^2/4 * J1(k*D/2) / (k*D/4)  = A_disk * 2*J1(kR)/(kR)
                    %    => m_hat(k) = pi*R^2 * 2*J1(kR)/(kR),   R = D/2
                    %
                    %    Then chi(r) = inverse Hankel of chi_hat(k)  [autocovariance]
                    %    S2(r) = phi^2 + chi(r)
                    % =========================================================
                    R     = D / 2;
                    kR    = k * R;
                    % Disk form factor in 2D (Fourier transform of indicator function of disk)
                    % m_hat(k) = 2*pi*R^2 * J1(kR)/(kR)  [area * normalized Bessel]
                    % Avoid division by zero at k=0
                    m_hat       = zeros(Nk, 1);
                    idx0        = kR < 1e-8;
                    m_hat(idx0) = pi * R^2;                        % limit as k->0
                    m_hat(~idx0)= 2*pi*R^2 * besselj(1, kR(~idx0)) ./ kR(~idx0);

                    % Spectral density (autocovariance in k-space)
                    chi_hat = rho * abs(m_hat).^2 .* S_k;

                    % Autocovariance chi(r) = S2(r) - phi^2  [inverse Hankel of chi_hat]
                    chi_r   = hankel_inv(chi_hat);

                    % Two-point correlation function
                    S2_r    = eta^2 + chi_r;

                    % Enforce physical bounds
                    S2_r    = max(0, min(1, S2_r));
                    Racf = (S2_r-eta^2)./(eta*(1-eta));
                    Racf = Racf - Racf(end);
                    Racf(2:end+1) = Racf;
                    Racf(1) = 1;
                    r(2:end+1) = r;
                    r(1) = 0;
                    obj.r = r;
                    obj.R = @(z) interp1(r, Racf, z, 'makima', 0);
                    obj.CalcLc;
                    r_norm = r / obj.CorrelationLength;
                    obj.r  = r_norm;
                    obj.R  = @(z) interp1(r_norm, Racf, z, 'makima', 0);
                    obj.CalcPhi;

                    out = obj.Phi;

                case 3
                   
                    % discretization in Fourier space
                    k = 0:0.1:3000;
                    r = 0:D/20:20*D;
                    % number density of spheres
                    rho = 3*eta/(4*pi);
                    % normalization of the radii
                    r = 2*r/D;
                    % Union volume of two spheres of unit radius
                    V2 = 8*pi/3*ones( size(r) );
                    V2( r<2 ) = 4*pi/3*(1+3/4*r(r<2)-r(r<2).^3/16);

                    % Fourier transform of the direct correlation function
                    % using the Percus-Yevick approximation
    
                    l1 = (1+2*eta)^2/(1-eta)^4;
                    l2 = -(1+eta/2)^2/(1-eta)^4;

                    c = -4*pi./(k.^3) .* ( l1*(sin(2*k)-2*k.*cos(2*k)) + ...
                        3*eta*l2./k.*( 4*k.*sin(2*k) + (2-4*k.^2).* ...
                        cos(2*k) - 2 ) + ...
                        eta*l1./(2*k.^3) .* (( 6*k.^2 - 3 -2*k.^4 ) .* ...
                        cos(2*k) + ...
                        (4*k.^3-6*k).*sin(2*k) + 3 ));
                    c(1) = -8*pi/3*((4+eta)*l1+18*eta*l2);

                    % Fourier transform of the total correlation function
                    % using the Ornstein-Zernike relation
                    h = c./(1-rho*c);
                    % h = c1./(1-rho*c1);
                    % Fourier transform of the Heaviside function

                    m = 4*pi * ( ((sin(k)./k) - cos(k))./(k.^2) );
                    m(1) = 4*pi/3;
                    % computation of 2-point matrix probability function
                    M = zeros( size(r) );
                    M(1) = -(eta/rho)^2;
                    for i1 = 2:length(r)
                        M(i1) = 1/2/pi^2/r(i1)*trapz( k, h.*m.^2.*k.*sin(k*r(i1)) );
                    end
                    S = 1 - rho*V2 + rho^2*M + eta^2;
                    % computation of the particle autocorrelation function
                    Racf = (S -(1-eta)^2) / eta / (1-eta);

                    % r is in D/2-normalised units; convert back to physical
                    % units so CalcLc produces a physical Lc, consistent with
                    % the 2D solver and with GetPSDFromImage output.
                    r_phys = r * D/2;
                    obj.r = r_phys;
                    obj.R = @(z) interp1(r_phys, Racf, z, 'makima', 0);
                    obj.CalcLc;
                    r_norm = r_phys / obj.CorrelationLength;
                    obj.r = [];
                    obj.R = [];
                    obj.r  = r_norm;
                    obj.R  = @(z) interp1(r_norm, Racf, z, 'makima', 0);
                    obj.CalcPhi;

                    out = obj.Phi;

            end
        end
        %% function to evaluate PSDF: an image
        function out = GetPSDFromImage(obj, Im, dx)
            %% GetPSDFromImage
            % Compute the normalised power spectral density function from a
            % 2D grayscale image or a 3D binary volume (logical/uint8/double).
            %
            %   1. Convert input to a binary (0/1) indicator field.
            %   2. Compute S2(r) via the Wiener-Khinchin theorem:
            %         S2_map = IFFT( |FFT(I)|^2 ) / N_voxels
            %   3. Radially average S2_map to obtain the isotropic S2(r).
            %   4. Derive the normalised autocorrelation (ACF):
            %         R(r) = [S2(r) - phi^2] / [phi*(1-phi)]
            %   5. Compute the correlation length Lc from R(r) in physical units.
            %   6. Normalise r by Lc and compute the PSDF via Fourier-Bessel
            %      quadrature. Store obj.Phi, obj.R, obj.k, obj.CorrelationLength.
            %
            % The resulting obj.CorrelationLength is in the same physical units
            % as dx, and is directly comparable with MonoDisperseSphere output
            % (which also stores Lc in physical units after the 3D PY fix).
            %
            % Syntax:
            %   out = GetPSDFromImage(obj, Im)
            %   out = GetPSDFromImage(obj, Im, dx)
            %
            % Inputs:
            %   Im : 2D or 3D numeric array (binary volume, values > 0.5 = solid).
            %   dx : isotropic voxel/pixel size in physical units (default 1)
            %
            % Output:
            %   out : function handle  Phi(k)  — normalised PSDF
            %
            % Side effects (properties set on obj):
            %   obj.Phi              — PSDF function handle
            %   obj.R                — normalised ACF function handle  R(r/Lc)
            %   obj.k                — wavenumber vector  k*Lc  [-]
            %   obj.r                — radial distance vector  r/Lc  [-]
            %   obj.CorrelationLength— physical correlation length  [same units as dx]
            %
            % See also: MaterialClass.CalcS2Correlation, MaterialClass.VoxelizeDomain

            % --- default pixel size ---
            if nargin < 3 || isempty(dx), dx = 1; end
            dy = dx;
            dz = dx;

            % --- convert to binary double indicator field ---
            if islogical(Im)
                I = single(Im);
            elseif ndims(Im) == 3 && size(Im,3) == 3
                % RGB image — convert to grayscale first
                I = single(rgb2gray(Im)) / 255;
                I = single(I > 0.5);
            else
                I = single(Im);
                if max(I(:)) > 1,  I = I / max(I(:));  end
                I = single(I > 0.5);
            end

            % =========================================================
            % 1D branch: vector input
            % =========================================================
            if isvector(I)
              error('Not implemented')
            end

            % =========================================================
            % 2D / 3D branch
            % =========================================================
            [nx, ny, nz] = size(I);
            if nz == 1
                % 2D image — set L so that min(L) = min physical dimension,
                % not a single voxel (which would limit max_r to half a pixel)
                L         = [nx*dx, ny*dy, min(nx*dx, ny*dy)];
                resolucao = min([dx dy]);
            else
                L         = [nx*dx, ny*dy, nz*dz];
                resolucao = min([dx dy dz]);
            end

            % --- S2 via FFT (CalcS2Correlation) ---
            [r_axis, S2_radial, ~, phi_vol, ~] = ...
                MaterialClass.CalcS2Correlation(I, resolucao, L);

            % --- normalised autocorrelation R(r) = [S2(r)-phi^2]/[phi*(1-phi)] ---
            denom = phi_vol * (1 - phi_vol);
            if denom < 1e-12
                warning('MaterialClass:GetPSDFromImage', ...
                    'Volume fraction is %.4f — image may be degenerate.', phi_vol);
                denom = 1;
            end
            if 1
                Racf   = (S2_radial - phi_vol^2) / denom;
                Racf   = Racf(:);
                Racf(2:end+1) = Racf;
                Racf(1) = 1;
                r_phys = r_axis(:);
                r_phys(2:end+1)=r_phys;
                r_phys(1) = 0;
            end
            if 0
                r_phys = r_axis(:);
                r_phys(2:end+1)=r_phys;
                r_phys(1) = 0;

                Racf = S2_radial - S2_radial(end);
                Racf = Racf./Racf(1);
                Racf(end+1) = 0;
            end
            if 0
                % correlation length — set obj.r/obj.R in physical units first so
            % CalcLc can integrate using the correct dimension
            obj.r = r_phys;
            obj.R = @(z) interp1(r_phys, Racf, z, 'makima', 0);
            Lc = obj.CalcLc;
            if Lc <= 0
                warning('MaterialClass:GetPSDFromImage', ...
                    'Correlation length came out non-positive (%.4g). Check image.', Lc);
                Lc = r_phys(end) / 10;
            end
            r_norm = r_phys / Lc;
            if any(r_norm > 20)
                idx = find(r_norm>20,1,'first');
                win = hann(21);
                win = [ones(idx,1); win(11:end); zeros(numel(r_norm)-idx-11,1)];
                Racf = Racf.*win;
            end
            end
            % Taper Racf: half-cosine from 1→0 between 80 % and 90 %, then 0.
            % 60-75% was too aggressive: for hard-sphere ACFs with correlation
            % length ~3D and max_r~15D, cutting at 9D lost ~40% of the
            % integral r^2*R(r), causing Lc to be underestimated by ~15-20%.
            N_racf  = numel(Racf);
            i_start = round(0.80 * N_racf);
            i_end   = round(0.90 * N_racf);
            n_tap   = i_end - i_start + 1;
            taper   = 0.5 * (1 + cos(pi * (0:n_tap-1)' / (n_tap - 1)));
            Racf(i_start:i_end) = Racf(i_start:i_end) .* taper;
            Racf(i_end+1:end)   = 0;

            obj.r = r_phys;
            obj.R = @(z) interp1(r_phys, Racf, z, 'makima', 0);
            Lc = obj.CalcLc;
            if ~isreal(Lc) || Lc <= 0
                figure
                plot(r_phys, Racf)
                warning('MaterialClass:GetPSDFromImage', ...
                    'Correlation length non-positive or complex (%.4g). Check image.', real(Lc));
                Lc = r_phys(end) / 10;
            end

            r_norm = r_phys / Lc;

            obj.r = [];
            obj.R = [];
            obj.r  = r_norm;
            obj.R  = @(z) interp1(r_norm, Racf, z, 'makima', 0);

            obj.CalcPhi;
            psd_vals = obj.Phi(obj.k);
            out     = obj.Phi;

            r_plot  = r_norm;
            S2_plot = S2_radial;
            S2_plot(end+1) = 0;
            R_plot  = Racf;
            dim_label = sprintf('%dD', 2 + (nz > 1));
            MaterialClass.psd_summary_figure(dim_label, phi_vol, Lc, ...
                r_plot, S2_plot, R_plot, obj.k, psd_vals, phi_vol);
        end
        function out = ImportPSDF(obj, Cin, rin)
            %% ImportPSDF
            % Import an externally computed autocorrelation function and
            % derive the normalised PSDF from it.
            %
            % Syntax:
            %   out = ImportPSDF(obj, Cin, rin)
            %
            % Inputs:
            %   Cin : autocorrelation values R(r),  vector of length N.
            %         Should start at R(0) ~ 1 and decay to 0.
            %   rin : corresponding radial distances r,  vector of length N.
            %         Must be in physical units (same as the desired Lc output).
            %         Must start at or near 0 and be monotonically increasing.
            %
            % Output:
            %   out : function handle  Phi(k)  — normalised 3D PSDF
            %
            % Side effects (properties set on obj):
            %   obj.Phi              — PSDF function handle
            %   obj.R                — normalised ACF function handle
            %   obj.r                — radial axis normalised by Lc
            %   obj.CorrelationLength— 3D correlation length Lc

            Cin = Cin(:);
            rin = rin(:);

            % --- 3D correlation length ---
            % Definition: lc^3 = 3 * int_0^inf r^2 * R(r) dr
            I3D  = 3 * trapz(rin, rin.^2 .* Cin);
            Lc3D = I3D^(1/3);
            if ~isreal(Lc3D) || Lc3D <= 0
                % Fallback to 1D definition: Lc = 2 * int_0^inf R(r) dr
                Lc3D = 2 * trapz(rin, Cin);
                warning('MaterialClass:ImportPSDF', ...
                    '3D Lc was non-positive; using 1D definition instead (Lc=%.4g).', Lc3D);
            end

            obj.CorrelationLength = Lc3D;
            obj.SpectralLaw       = 'Imported';

            % --- normalise r by Lc ---
            rlc   = rin / Lc3D;
            obj.r = rlc;
            obj.R = @(z) interp1(rlc, Cin, z, 'pchip', 0);

            % --- 3D PSDF via Fourier-Bessel (sinc) quadrature ---
            %   Phi(k) = 4*pi/(2*pi)^3 * int_0^inf r^2 sinc(k*r) R(r) dr
            %   with r and k normalised by Lc
            k = logspace(log10(1e-3), log10(6), 512);
            PSD = zeros(size(k));
            for i = 1:length(k)
                if k(i) < 1e-8
                    integrand = rlc.^2 .* Cin;
                else
                    kr = k(i) * rlc;
                    j0 = sin(kr) ./ kr;
                    j0(kr < 1e-12) = 1;
                    integrand = rlc.^2 .* Cin .* j0;
                end
                PSD(i) = (4*pi / (2*pi)^3) * trapz(rlc, integrand);
            end
            PSD = max(PSD, 0);    % enforce non-negativity

            obj.k   = k;
            obj.Phi = @(z) interp1(k, PSD, z, 'pchip', 0);
            out     = obj.Phi;
        end
        %% Calc
        function Lc = CalcLc(obj)
            %% CalcLc
            % Compute the correlation length from the correlation function
            %
            % Syntax:
            %   Lc = CalcLc (  );
            %
            % Inputs:
            %   (none - uses obj.R function handle)
            %
            % Output:
            %   Lc: correlation length scalar
            %
            % The correlation functions are normalized in 1D, 2D, 3D as follows:
            %   1D : lc = 2 * integral_0^inf R(x) dx
            %   2D : lc^2 = 2 * integral_0^inf x*R(x) dx
            %   3D : lc^3 = 3 * integral_0^inf x^2*R(x) dx
            %
            % For isotropic correlation functions in d dimensions.
            
            if isempty(obj.r)
                r = linspace(0, 15, 4096);
                r(1) = 1e-6;
            else
                r = obj.r;
            end
            R_vals = obj.R(r);

            if obj.d == 1
                Lc = abs(2 * trapz(r, R_vals));
            elseif obj.d == 2
                I  = trapz(r, r .* R_vals);
                Lc = sqrt(2 * abs(I));
            elseif obj.d == 3
                I  = trapz(r, r.^2 .* R_vals);
                Lc = nthroot(3 * abs(I), 3);
            else
                error('MaterialClass:InvalidDimension', 'Dimension must be 1, 2, or 3');
            end
            obj.CorrelationLength = Lc;
        end
        function out = CalcPhi(obj)
            %% CalcPhi
            % Compute the power spectral density from the correlation function
            %
            % Syntax:
            %   out = CalcPhi (  );
            %
            % Inputs:
            %   (none - uses obj.R function handle)
            %
            % Output:
            %   out: function handle Phi(k) for the power spectral density
            %
            % The power spectral density is the Fourier transform of the 
            % correlation function. For isotropic functions in d dimensions:
            %   1D : Phi(k) = 2/(2*pi) * integral_0^inf cos(k*x)*R(x) dx
            %   2D : Phi(k) = 2*pi/(2*pi)^2 * integral_0^inf x*J_0(k*x)*R(x) dx
            %   3D : Phi(k) = 4*pi/(2*pi)^3 * integral_0^inf x^2*sinc(k*x)*R(x) dx
            %
            % where sinc(k*x) = sin(k*x)/(k*x)
            
            % Grid for real space (r) and wavenumber (k)
            r = linspace(0, 15, 4096);  % real space grid (normalised by Lc)
            r(1) = 1e-6;                % avoid r=0 in sinc/Bessel kernels
            R_vals = obj.R(r);

            % Wavenumber grid
            k = linspace(0, 20, 4096);  % wavenumber grid (normalised by Lc)
            k(1) = 1e-6;                % avoid k=0 in sinc kernel
            PSD = zeros(size(k));
            
            if obj.d == 1
                % Phi(k) = 2/(2*pi) * int_0^inf cos(k*r)*R(r) dr
                for i = 1:length(k)
                    integrand = cos(k(i) * r) .* R_vals;
                    PSD(i) = (2 / (2*pi)) * trapz(r, integrand);
                end
            elseif obj.d == 2
                % Phi(k) = 2*pi/(2*pi)^2 * int_0^inf r*J_0(k*r)*R(r) dr
                for i = 1:length(k)
                    J0 = besselj(0, k(i) * r);
                    integrand = r .* J0 .* R_vals;
                    PSD(i) = (2*pi / (2*pi)^2) * trapz(r, integrand);
                end
            elseif obj.d == 3
                % Phi(k) = 4*pi/(2*pi)^3 * int_0^inf r^2*sinc(k*r)*R(r) dr
                % where sinc(k*r) = sin(k*r)/(k*r)
                for i = 1:length(k)
                    kr = k(i) * r;
                    j0 = sin(kr) ./ kr;
                    j0(kr < 1e-12) = 1;  % handle sinc(0) = 1
                    integrand = r.^2 .* j0 .* R_vals;
                    PSD(i) = (4*pi / (2*pi)^3) * trapz(r, integrand);
                end
            else
                error('MaterialClass:InvalidDimension', 'Dimension must be 1, 2, or 3');
            end
            
            % Ensure non-negativity (PSD should be >= 0)
            PSD = max(PSD, 0);
            
            % Store results in object
            %obj.r = r;
            obj.k = k;
            obj.Phi = @(z) interp1(k, PSD, z, 'makima', 0);
            out = obj.Phi;
        end
        function out = CalcR(obj)
            %% CalcR
            % Compute the correlation function from the power spectral density
            %
            % Syntax:
            %   out = CalcR (  );
            %
            % Inputs:
            %   (none - uses obj.Phi function handle)
            %
            % Output:
            %   out: function handle R(x) for the correlation function
            %
            % This is the inverse Fourier transform of the PSD.
            % For isotropic functions in d dimensions:
            %   1D : R(x) = 2*pi * integral_0^inf cos(k*x)*Phi(k) dk
            %   2D : R(x) = 2*pi * integral_0^inf k*J_0(k*x)*Phi(k) dk
            %   3D : R(x) = (2*pi)^3/(4*pi) * integral_0^inf k^2*sinc(k*x)*Phi(k) dk
            %
            % where sinc(k*x) = sin(k*x)/(k*x)
            
            % Grid for wavenumber (k) and real space (r)
            k = logspace(-3, 2, 512);  % wavenumber grid
            Phi_vals = obj.Phi(k);
            
            % Real space grid
            r = logspace(-4, 3, 4096);  % real space grid
            R_vals = zeros(size(r));
            
            if obj.d == 1
                % R(r) = 2*pi * int_0^inf cos(k*r)*Phi(k) dk
                for i = 1:length(r)
                    integrand = cos(k * r(i)) .* Phi_vals;
                    R_vals(i) = 2*pi * trapz(k, integrand);
                end
            elseif obj.d == 2
                % R(r) = 2*pi * int_0^inf k*J_0(k*r)*Phi(k) dk
                for i = 1:length(r)
                    J0 = besselj(0, k * r(i));
                    integrand = k .* J0 .* Phi_vals;
                    R_vals(i) = 2*pi * trapz(k, integrand);
                end
            elseif obj.d == 3
                % R(r) = (2*pi)^3/(4*pi) * int_0^inf k^2*sinc(k*r)*Phi(k) dk
                % where sinc(k*r) = sin(k*r)/(k*r)
                for i = 1:length(r)
                    kr = k * r(i);
                    j0 = sin(kr) ./ kr;
                    j0(kr < 1e-12) = 1;  % handle sinc(0) = 1
                    integrand = k.^2 .* j0 .* Phi_vals;
                    R_vals(i) = ((2*pi)^3 / (4*pi)) * trapz(k, integrand);
                end
            else
                error('MaterialClass:InvalidDimension', 'Dimension must be 1, 2, or 3');
            end
            
            % Normalize R(0) = 1 (correlation function should start at 1)
            R_vals = R_vals / R_vals(1);
            
            % Store results in object
            obj.k = k;
            obj.r = r;
            obj.R = @(z) interp1(r, R_vals, z, 'pchip', 0);
            out = obj.R;
        end
        %% plot
        function h = PlotPSD(obj,h)
            if ~exist('h','var')
                h = figure;
            end
            if isempty(obj.k)
                obj.k = linspace(0,10,1024);
            end
            plot(obj.k,obj.Phi(obj.k),'LineWidth',2)
            xlabel('Normalized wavenumer [-]')
            ylabel('Power Spectral Density [-]')
            saux = sprintf("PSD model:%s, dimension:%d",obj.SpectralLaw,obj.d);
            title(saux)
            grid on
            box on
            set(gca,'FontSize',14)
        end
        function h = PlotCorrelation(obj,h)
            if ~exist('h','var')
                h = figure;
            end
            x = linspace(0,10,2048);
            plot(x,obj.R(x),'LineWidth',2)
            xlabel('Normalized lag distance [-]')
            ylabel('Correlation [-]')
            grid on
            box on
            saux = sprintf("Correlation model:%s, dimension:%d",obj.SpectralLaw,obj.d);
            title(saux)
            set(gca,'FontSize',14)
        end
        function h = plotsigma(obj,h)
            if ~exist('h','var')
                h = figure;
            end
            z = linspace(0,2*pi,2*2048);
            if obj.acoustics
                hold on
                plot(z,obj.sigma{1}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('Differential Scattering Cross-Sections [-]')
                grid on
                box on
                saux = sprintf("Model:%s, dimension:%d",obj.SpectralLaw,obj.d);
                title(saux)
                set(gca,'FontSize',14)
            else
                subplot(2,2,1)
                hold on
                plot(z,obj.sigma{1,1}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('P2P [-]')
                grid on
                box on
                set(gca,'FontSize',14)

                subplot(2,2,2)
                hold on
                plot(z,obj.sigma{1,2}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('P2S [-]')
                grid on
                box on
                set(gca,'FontSize',14)

                subplot(2,2,3)
                hold on
                plot(z,obj.sigma{2,1}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('S2P [-]')
                grid on
                box on
                set(gca,'FontSize',14)

                subplot(2,2,4)
                hold on
                plot(z,obj.sigma{2,2}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('S2S [-]')
                grid on
                box on
                set(gca,'FontSize',14)
                saux = sprintf("Model:%s, dimension:%d",obj.SpectralLaw,obj.d);
                sgtitle(saux)
            end
        end
        function h = plotpolarsigma(obj,h)
            if ~exist('h','var')
                h = figure;
            end
            z = linspace(0,2*pi,2*2048);
            if obj.acoustics
                polaraxes
                hold on
                polarplot(z,obj.sigma{1}(z),'LineWidth',2)
                %xlabel('Angle [rad]')
                title('Differential Scattering Cross-Sections')
                grid on
                box on
                set(gca,'FontSize',14)
            else
                subplot(2,2,1)
                polarplot(z,obj.sigma{1,1}(z),'LineWidth',2)
                title('P2P [-]')
                grid on
                box on
                set(gca,'FontSize',14)

                subplot(2,2,2)
                polarplot(z,obj.sigma{1,2}(z),'LineWidth',2)
                title('P2S [-]')
                grid on
                box on
                set(gca,'FontSize',14)

                subplot(2,2,3)
                polarplot(z,obj.sigma{2,1}(z),'LineWidth',2)
                title('S2P [-]')
                grid on
                box on
                set(gca,'FontSize',14)

                subplot(2,2,4)
                polarplot(z,obj.sigma{2,2}(z),'LineWidth',2)
                title('S2S [-]')
                grid on
                box on
                set(gca,'FontSize',14)
            end
        end
        %% TRANSMISSION REFLECTION
        function o = MaterialInterface(obj, P, n, ind)
            % MaterialInterface
            % Compute reflection (and optionally transmission) of particles
            % at an interface. This implementation assumes a solid–air
            % interface, thus pure reflection and no transmission.
            %
            % P   - structure of particles (fields: p, dir, perp, ...)
            % n   - interface normal (1x3)
            % ind - particles outside the domain (those that hit the boundary)

            if isempty(obj.rho)
                error(['Please define the density in MaterialClass before ' ...
                    'calling MaterialInterface.']);
            end

            % Split P vs S particles at the boundary
            pIdx =  P.p  & ind;   % incident P-waves
            sIdx = ~P.p  & ind;   % incident S-waves

            % --- P-WAVES PROCESSING ---
            partP_dir  = P.dir(pIdx,:);
            nP = size(partP_dir,1);

            if nP > 0
                nRepP = repmat(n, nP, 1);
                dotPn = dot(partP_dir, nRepP, 2);
                absDotPn = abs(dotPn);
                absDotPn(absDotPn > 1) = 1;
                incAngP = acosd(absDotPn);
            else
                incAngP = [];
            end

            % --- S-WAVES PROCESSING (SV/SH SPLIT) ---
            partS_dir  = P.dir(sIdx,:);
            partS_perp = P.perp(sIdx,:);
            nS = size(partS_dir,1);

            if nS > 0
                nRepS = repmat(n, nS, 1);
                % SH direction: orthogonal to plane of incidence {k, n}
                eSH = cross(partS_dir, nRepS, 2);
                norm_eSH = vecnorm(eSH, 2, 2);
                % Handle near-normal incidence (k || n)
                nearNormal = norm_eSH < eps;

                if any(nearNormal)
                    % Select the coordinate axis least aligned with the normal.
                    nUnit = n(:).';
                    [~,referenceIndex] = min(abs(nUnit));
                    referenceAxis = zeros(size(nUnit));
                    referenceAxis(referenceIndex) = 1;
                    t = repmat(referenceAxis,sum(nearNormal),1);
                    t = t - (t * nUnit.') * nUnit;
                    eSH(nearNormal,:)    = t;
                    norm_eSH(nearNormal) = vecnorm(t, 2, 2);
                end

                eSH = eSH ./ norm_eSH;

                % SV direction: in plane {k,n}, orthogonal to k
                eSV = cross(eSH, partS_dir, 2);
                eSV = eSV ./ vecnorm(eSV, 2, 2);

                % Project polarization onto SV/SH basis
                aSV = dot(partS_perp, eSV, 2);
                aSH = dot(partS_perp, eSH, 2);

                % Energy fractions
                E_SV = aSV.^2;
                E_SH = aSH.^2;
                E_tot = E_SV + E_SH + eps;

                % Probability of being SV
                probSV = E_SV ./ E_tot;

                % Monte Carlo decision: SV or SH?
                isSV = rand(size(probSV)) < probSV;

                % Calculate Incidence Angle for S
                dotSn   = dot(partS_dir, nRepS, 2);
                absDotSn = abs(dotSn);
                absDotSn(absDotSn > 1) = 1; % Clamp

                incAngS = acosd(absDotSn);

                incAngSV = incAngS(isSV);
                incAngSH = incAngS(~isSV);

                % Pointers to global indices
                partSV_dir = partS_dir(isSV, :);
                partSH_dir = partS_dir(~isSV, :);
            else
                isSV = false(0,1);
                incAngSV = [];
                incAngSH = [];
                partSV_dir = [];
                partSH_dir = [];
            end

            % --- REFLECTION COEFFICIENTS ---
            [out, angles] = obj.getZoeppritzCached();

            % Helper for fast lookup
            lookup = @(data, ang) interp1(out.j1_deg, data, ang, 'linear', 'extrap');

            % --- P INCIDENCE ---
            if nP > 0
                % Get Energy Coefficient for P -> P reflection
                Rpp_coeff = lookup(out.E_Rpp, incAngP);

                p2p  = rand(nP,1) <= Rpp_coeff; % P -> P reflection
                p2sv = ~p2p;                    % P -> SV reflection

                % Get Angles (Snell's law is embedded in Zoeppritz output or recomputed)
                angP2P  = angles.Rpp(incAngP(p2p))';
                angP2SV = angles.Rpsv(incAngP(p2sv))';

                % Get Amplitudes
                ampP2P  = lookup(out.A_Rpp, incAngP(p2p));
                ampP2SV = lookup(out.A_Rpsv, incAngP(p2sv));
            else
                p2p = []; p2sv = [];
                angP2P = []; angP2SV = [];
                ampP2P = []; ampP2SV = [];
            end

            % --- SV INCIDENCE ---
            nSV = size(incAngSV, 1);
            if nSV > 0
                Rsvsv_coeff = lookup(out.E_Rsvsv, incAngSV);

                sv2sv = rand(nSV,1) <= Rsvsv_coeff;
                sv2p  = ~sv2sv;

                angSV2SV = angles.Rsvsv(incAngSV(sv2sv))';
                angSV2P  = angles.Rsvp(incAngSV(sv2p))';

                ampSV2SV = lookup(out.A_Rsvsv, incAngSV(sv2sv));
                ampSV2P  = lookup(out.A_Rsp, incAngSV(sv2p));
            else
                sv2sv = []; sv2p = [];
                angSV2SV = []; angSV2P = [];
                ampSV2SV = []; ampSV2P = [];
            end

            % --- SH INCIDENCE ---
            nSH = size(incAngSH, 1);
            if nSH > 0
                % Free surface: SH reflects as SH (R=1)
                sh2sh = true(nSH, 1);

                angSH2SH = angles.Rsh(incAngSH)';
                ampSH2SH = lookup(out.A_Rsh, incAngSH);
            else
                sh2sh = [];
                angSH2SH = [];
                ampSH2SH = [];
            end

            % --- OUTPUT ---
            o.pparticle  = pIdx;
            o.sparticle  = sIdx;
            o.svparticle = isSV;
            o.shparticle = ~isSV;

            o.p2p   = p2p;
            o.p2sv  = p2sv;

            o.sv2sv = sv2sv;
            o.sv2p  = sv2p;

            o.sh2sh = sh2sh;

            % Scattering angles (optional, if solver calculates them itself)
            o.angPP   = angP2P;
            o.angPSV  = angP2SV;
            o.angSVSV = angSV2SV;
            o.angSVP  = angSV2P;
            o.angSHSH = angSH2SH;

            o.ampPP   = ampP2P;
            o.ampPSV  = ampP2SV;
            o.ampSVSV = ampSV2SV;
            o.ampSVP  = ampSV2P;
            o.ampSHSH = ampSH2SH;
        end

        function [out, angles] = getZoeppritzCached(obj)
            % Reuse the coefficients while the material properties remain unchanged.
            cacheKey = [obj.vp obj.vs obj.rho obj.acoustics];

            if isempty(obj.zoeppritzOutputCache) || isempty(obj.zoeppritzCacheKey) || ...
                    ~isequaln(obj.zoeppritzCacheKey,cacheKey)
                [obj.zoeppritzOutputCache,obj.zoeppritzAnglesCache] = MaterialClass.Zoeppritz(obj);
                obj.zoeppritzCacheKey = cacheKey;
            end

            out = obj.zoeppritzOutputCache;
            angles = obj.zoeppritzAnglesCache;
        end
    end
    methods(Static)
        function material = polycrystal3D(geometry,frequency,crystalInput, ...
                tpcfInput,tpcfParameter,phaseVelocities)
            %% polycrystal3D
            % Create a 3D elastic MaterialClass object for a polycrystal.
            %
            % Monodisperse exponential or spherical TPCF:
            %   material = MaterialClass.polycrystal3D(geometry,frequency, ...
            %       crystalInput,'exp',correlationLength)
            %   material = MaterialClass.polycrystal3D(geometry,frequency, ...
            %       crystalInput,'spherical',grainDiameter)
            %
            % Grain size distributions (lognormal, Gamma, Weibull, or
            % normal truncated at zero):
            %   material = MaterialClass.polycrystal3D(geometry,frequency, ...
            %       crystalInput,'ShengKhazaie',grainSizeDistribution)
            %
            % User-defined dimensional spectral TPCF etaTilde(q):
            %   material = MaterialClass.polycrystal3D(geometry,frequency, ...
            %       crystalInput,etaTilde)
            %
            % In either form, an optional final [Vp Vs] input overrides the
            % Voigt velocities calculated from the single-crystal stiffness.

            if nargin < 6, phaseVelocities = []; end
            if nargin < 5, tpcfParameter = []; end

            if ~isstruct(geometry) || ~isfield(geometry,'dimension') || ...
                    ~isscalar(geometry.dimension) || geometry.dimension ~= 3
                error('MaterialClass:PolycrystalGeometry', ...
                      'MaterialClass.polycrystal3D requires a 3D geometry.');
            end
            validateattributes(frequency, {'numeric'}, ...
                {'real','finite','scalar','positive'}, ...
                 'MaterialClass.polycrystal3D', 'frequency');

            if ischar(crystalInput) || (isstring(crystalInput) && isscalar(crystalInput))
                crystal = MaterialClass.singleCrystalProperties(crystalInput);
            elseif isstruct(crystalInput) && isscalar(crystalInput) && ...
                    isfield(crystalInput,'rho') && isfield(crystalInput,'C')
                crystal = crystalInput;
                if ~isfield(crystal,'name'), crystal.name = 'Custom'; end
                if ~isfield(crystal,'symmetry')
                    crystal.symmetry = 'unspecified';
                end
                if ~isfield(crystal,'constants')
                    crystal.constants = struct();
                end
            else
                error('MaterialClass:InvalidSingleCrystalInput', ...
                    ['crystalInput must be a material name or a scalar ', ...
                     'structure containing rho and C.']);
            end

            if isa(tpcfInput,'function_handle')
                if nargin > 5
                    error('MaterialClass:TooManyCustomTPCFInputs', ...
                         ['With a user-defined spectral TPCF, append at most ', ...
                          'one optional [Vp Vs] input.']);
                end
                phaseVelocities = tpcfParameter;
                TPCF = struct('model','userDefined', ...
                              'parameters',struct(), ...
                              'spectrum',tpcfInput);
            elseif ischar(tpcfInput) || (isstring(tpcfInput) && isscalar(tpcfInput))
                if isempty(tpcfParameter)
                    error('MaterialClass:MissingTPCFParameter', ...
                          'The selected TPCF model requires a model parameter.');
                end
                TPCF = MaterialClass.PolycrystalTPCF(tpcfInput,tpcfParameter);
            else
                error('MaterialClass:InvalidPolycrystalTPCFInput', ...
                     ['tpcfInput must be a model name or a function handle ', ...
                      'for the dimensional spectral TPCF.']);
            end

            if ~isempty(phaseVelocities)
                if ~isnumeric(phaseVelocities) || ...
                        ~isreal(phaseVelocities) || ...
                        numel(phaseVelocities) ~= 2 || ...
                        any(~isfinite(phaseVelocities)) || ...
                        any(phaseVelocities <= 0)
                    error('MaterialClass:InvalidPolycrystalVelocities', ...
                        'phaseVelocities must be positive [Vp Vs].');
                end
                phaseVelocities = reshape(phaseVelocities,1,2);
            end

            material = MaterialClass();
            material.d = geometry.dimension;
            material.acoustics = false;
            material.Frequency = frequency;
            material.scatteringModel = 'polycrystal';
            material.correlationStructure = 'isotropic';
            material.singleCrystal = crystal;
            material.rho = crystal.rho;
            material.TPCF = TPCF;
            if ~isempty(phaseVelocities)
                material.vp = phaseVelocities(1);
                material.vs = phaseVelocities(2);
            end

            % Construct only the differential scattering cross-sections.
            % prepareSigma is called through the normal solver pathway.
            material.CalcSigma;
        end
        function material = singleCrystalProperties(materialName)
            %% singleCrystalProperties
            % Return named single-crystal density and elastic properties.
            %
            % Material names are case-insensitive. Spaces, hyphens, and
            % underscores are ignored. The returned structure contains the
            % canonical name, symmetry, density [kg/m^3], independent
            % constants [Pa], and 6-by-6 Voigt stiffness matrix C [Pa].
            % The property database was adapted from the material list
            % developed by Ningyue Sheng during his PhD research.

            if ~(ischar(materialName) || ...
                    (isstring(materialName) && isscalar(materialName)))
                error('MaterialClass:InvalidSingleCrystalName', ...
                    ['materialName must be a character vector or a ', ...
                     'scalar string.']);
            end

            key = lower(strtrim(char(materialName)));
            key = regexprep(key,'[\s_-]',''); % Remove spaces, underscores, and hyphens

            switch key
                % cubic material
                case {'al','aluminum','aluminium'}
                    material = MaterialClass.cubicSingleCrystal('Aluminum',2700,108,62,28.3);
                case {'cr','chromium'}
                    material = MaterialClass.cubicSingleCrystal('Chromium',7150,348,67,100);
                case {'nb','niobium'}
                    material = MaterialClass.cubicSingleCrystal('Niobium',8570,245,132,28.4);
                case {'pt','platinum'}
                    material = MaterialClass.cubicSingleCrystal('Platinum',21450,347,251,76.5);
                case {'fe','afe','alphafe','alphairon','ironalpha','ferrite'}
                    material = MaterialClass.cubicSingleCrystal('Alpha iron',7800,231,135,115);
                case {'la','lanthanum','lanthanium'}
                    material = MaterialClass.cubicSingleCrystal('Lanthanum',6145,34.5,20.4,18);
                case {'ni','nickel'}
                    material = MaterialClass.cubicSingleCrystal('Nickel',8900,247,153,122);
                case {'au','gold'}
                    material = MaterialClass.cubicSingleCrystal('Gold',19300,191,162,42.2);
                case {'ag','silver'}
                    material = MaterialClass.cubicSingleCrystal('Silver',10490,122,92,45.5);
                case {'cobaltcubic','cubiccobalt','fccco','cofcc'}
                    material = MaterialClass.cubicSingleCrystal('Cobalt (cubic)',8900,242,160,128);
                case {'cu','copper'}
                    material = MaterialClass.cubicSingleCrystal('Copper',8960,168.4,121.4,75.39);
                case {'pb','lead'}
                    material = MaterialClass.cubicSingleCrystal('Lead',11340,48.8,41.4,14.8);
                case {'gfe','gammafe','gammairon','irongamma','austenite'}
                    material = MaterialClass.cubicSingleCrystal('Gamma iron',8000,154,122,77);
                case {'k','potassium'}
                    material = MaterialClass.cubicSingleCrystal('Potassium',890,3.71,3.15,1.88);
                case {'li','lithium'}
                    material = MaterialClass.cubicSingleCrystal('Lithium',534,13.4,11.3,9.6);
                case {'inconel','inconnel'}
                    material = MaterialClass.cubicSingleCrystal('Inconel',8260,234.6,145.9,126.2);
                % hexagonal material
                case {'ti','ati','alphati','titanium','alphatitanium','titaniumalpha'}
                    material = MaterialClass.hexagonalSingleCrystal('Alpha titanium',4500,170,92,70,192,52);
                case {'zr','zirconium'}
                    material = MaterialClass.hexagonalSingleCrystal('Zirconium',6490,136,78,68,163,40);
                case {'be','beryllium'}
                    material = MaterialClass.hexagonalSingleCrystal('Beryllium',1850,292.3,26.7,14,336.4,162.5);
                case {'cd','cadmium'}
                    material = MaterialClass.hexagonalSingleCrystal('Cadmium',8650,114.1,41,40.3,49.9,19);
                case {'cobalthexagonal','hexagonalcobalt','cobalthex','hcpco','cohcp'}
                    material = MaterialClass.hexagonalSingleCrystal('Cobalt (hexagonal)',8900,295,159,111,335,71);
                case {'gd','gadolinium'}
                    material = MaterialClass.hexagonalSingleCrystal('Gadolinium',7900,66.7,25,21.3,71.9,20.7);
                case {'ho','holmium'}
                    material = MaterialClass.hexagonalSingleCrystal('Holmium',8800,76.5,25.6,21,79.6,25.9);
                case {'mg','magnesium'}
                    material = MaterialClass.hexagonalSingleCrystal('Magnesium',1738,59.3,25.7,21.4,61.5,16.4);
                case {'nd','neodymium'}
                    material = MaterialClass.hexagonalSingleCrystal('Neodymium',7010,54.8,24.6,16.6,60.9,15);
                case {'re','rhenium'}
                    material = MaterialClass.hexagonalSingleCrystal('Rhenium',21020,616,273,206,683,161);
                case {'ru','ruthenium'}
                    material = MaterialClass.hexagonalSingleCrystal('Ruthenium',12200,563,188,168,624,181);
                case {'sc','scandium'}
                    material = MaterialClass.hexagonalSingleCrystal('Scandium',3000,99.3,39.7,29.4,107,27.7);
                case {'tb','terbium'}
                    material = MaterialClass.hexagonalSingleCrystal('Terbium',8300,67.9,24.3,23,72.2,21.4);
                case {'tl','thallium'}
                    material = MaterialClass.hexagonalSingleCrystal('Thallium',11710,40.8,35.4,29,52.8,7.26);
                case {'y','yttrium'}
                    material = MaterialClass.hexagonalSingleCrystal('Yttrium',4472,77.9,29.2,20,76.9,24.3);
                case {'zn','zinc'}
                    material = MaterialClass.hexagonalSingleCrystal('Zinc',7133,165,31.1,50,61.8,39.6);
                case {'co','cobalt'}
                    error('MaterialClass:AmbiguousCobalt', ...
                         ['Cobalt is available in cubic and hexagonal forms. ', ...
                          'Use ''cobaltCubic'' or ''cobaltHexagonal''.']);
                otherwise
                    error('MaterialClass:UnknownSingleCrystal', ...
                         ['Unknown material "%s". Examples are ''Al'', ', ...
                          '''alphaIron'', ''Cu'', ''Ni'', ''Ti'', and ''Zr''.'], ...
                        char(materialName));
            end
        end
    end
    methods(Static, Access=private)
        [L,M,N] = InnerProducts(Cloc)
        TPCF = PolycrystalTPCF(modelName,modelParameter)
        function material = cubicSingleCrystal(name,rho,c11,c12,c44)
            % Construct a cubic single-crystal property structure
            scale = 1e9;
            c11 = c11*scale;
            c12 = c12*scale;
            c44 = c44*scale;
            C = [c11 c12 c12 0   0   0; ...
                 c12 c11 c12 0   0   0; ...
                 c12 c12 c11 0   0   0; ...
                 0   0   0   c44 0   0; ...
                 0   0   0   0   c44 0; ...
                 0   0   0   0   0   c44];
            constants = struct('c11',c11,'c12',c12,'c44',c44);
            material = struct('name',name,'symmetry','cubic', ...
                'rho',rho,'constants',constants,'C',C);
        end
        function material = hexagonalSingleCrystal( ...
                name,rho,c11,c12,c13,c33,c44)
            % Construct a hexagonal single-crystal property structure
            scale = 1e9;
            c11 = c11*scale;
            c12 = c12*scale;
            c13 = c13*scale;
            c33 = c33*scale;
            c44 = c44*scale;
            c66 = 0.5*(c11-c12);
            C = [c11 c12 c13 0   0   0; ...
                 c12 c11 c13 0   0   0; ...
                 c13 c13 c33 0   0   0; ...
                 0   0   0   c44 0   0; ...
                 0   0   0   0   c44 0; ...
                 0   0   0   0   0   c66];
            constants = struct('c11',c11,'c12',c12,'c13',c13,'c33',c33,'c44',c44,'c66',c66);
            material = struct('name',name,'symmetry','hexagonal','rho',rho,'constants',constants,'C',C);
        end
    end
    methods(Static)
        %% OTHER
        function h = heaviside(h)
            h(h>0) = 1;
            h(h==0) = 1/2;
            h(h<0) = 0;
        end
        function obj = preset(n)
            obj = MaterialClass();
            switch n
                case 1
                    obj.sigma{1} = @(th) 1/4/pi*ones(size(th));
                    obj.v = 2;
                    obj.acoustics = true;
                case 2
                    obj.sigma{1} = @(th) 1/10/pi*ones(size(th));
                    obj.v = 2;
                    obj.acoustics = true;
                case 3
                    obj.vp = 6;
                    obj.vs = 6/sqrt(3);
                    obj.acoustics = false;
                case 4
                    obj.vp = 6;
                    obj.vs = 6/sqrt(3);
                    obj.acoustics = false;
            end
        end
        function [out, angles] = Zoeppritz(mat)
            %% Zoeppritz
            % function to calculate reflection and transmission
            % coefficients for P, SV and SH waves
            %
            % Syntax:
            %   MaterialClass.Zoeppritz( mat );
            %
            % Inputs:
            %  mat: scalar elastic MaterialClass object. The exterior
            %       medium is currently hard-coded as air.
            %
            % Output: The method return the coefficients of transmission
            % and reflection or only reflection

            if isscalar(mat)
                % assuming solid-fluid interface
                % the fluid is hard coded as air
                vpAir  = 343; %m/s
                vsAir  = 0;   %m/s
                rhoAir = 1;  %kg/m^3

                % angles
                j1_deg = linspace(0,90,181);

                if(mat.acoustics)
                    error("This dont work for acoustic material")
                end

                vp1  = mat(1).vp;
                vs1  = mat(1).vs;
                rho1 = mat(1).rho;

                vp2  = vpAir;
                vs2  = vsAir;
                rho2 = rhoAir;

                if vs2 == 0
                    out = MaterialClass.ZoeppritzFluid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
                else
                    out = MaterialClass.ZoeppritzSolid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
                end

                out.j1_deg = j1_deg;

                % angles P
                angles.Rpp  = @(theta1) asind((sind(theta1) / vp1) * vp1);
                angles.Rpsv = @(theta1) asind((sind(theta1) / vp1) * vs1);
                angles.Tpp  = @(theta1) asind((sind(theta1) / vp1) * vp2);
                angles.Tpsv = @(theta1) asind((sind(theta1) / vp1) * vs2);

                % angles SV
                angles.Rsvsv = @(theta1) asind(sind(theta1) / vs1 * vs1);
                angles.Rsvp  = @(theta1) asind(sind(theta1) / vs1 * vp1);
                angles.Tsvsv = @(theta1) asind(sind(theta1) / vs1 * vs2);
                angles.Tsvp  = @(theta1) asind(sind(theta1) / vs1 * vp2);

                % Incident SH-wave
                angles.Rsh = @(theta1) asind(sind(theta1) / vs1 * vs1);
                angles.Tsh = @(theta1) asind(sind(theta1) / vs1 * vs2);

            else
                error('MaterialClass:MaterialInterfacesNotImplemented', ...
                    ['Zoeppritz currently accepts only one material and ', ...
                     'models its interface with air. Interfaces between ', ...
                     'multiple user-defined materials are not implemented.']);
            end
        end
        out = ZoeppritzFluid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
        out = ZoeppritzSolid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
        function Microstruture = VoxelizeDomain(coordinates, D, L, resolution)
            %% VoxelizeDomain
            % Discretise a periodic domain into a binary pixel/voxel grid and
            % paint rods (1D), disks (2D), or spheres (3D) at the given centres.
            %
            % Syntax:
            %   M = MaterialClass.VoxelizeDomain(coordinates, D, L, resolution)
            %
            % Inputs:
            %   coordinates : N x d  array of object centres  (d = 1, 2, or 3)
            %   D           : object diameter (scalar)
            %   L           : domain size — scalar (1D), 1x2 (2D), or 1x3 (3D)
            %   resolution  : pixel/voxel edge length (isotropic)
            %
            % Output:
            %   Microstruture : logical array
            %                   1D → nx x 1
            %                   2D → nx x ny
            %                   3D → nx x ny x nz
            %                   true = solid phase, false = matrix phase
            %
            % Examples:
            %   % 1D rods
            %   c1 = rand(50,1) * 1;
            %   M1 = MaterialClass.VoxelizeDomain(c1, 0.05, 1, 0.002);
            %
            %   % 2D disks
            %   c2 = rand(100,2) .* [1 1];
            %   M2 = MaterialClass.VoxelizeDomain(c2, 0.05, [1 1], 0.002);
            %
            %   % 3D spheres
            %   c3 = rand(200,3) .* [100 100 100];
            %   M3 = MaterialClass.VoxelizeDomain(c3, 4, [100 100 100], 1);

            dim = size(coordinates, 2);
            R   = D / 2;
            n_obj = size(coordinates, 1);
            report_step = max(1, floor(n_obj / 10));

            fprintf('Voxelizando o dominio (%dD)...\n', dim);

            switch dim
                % ---------------------------------------------------------
                case 1
                    nx = ceil(L / resolution);
                    xv = linspace(0, L, nx);
                    Microstruture = false(nx, 1);

                    for i = 1:n_obj
                        if n_obj > 100 && mod(i, report_step) == 0
                            fprintf('  Progresso: %.1f%%\n', 100*i/n_obj);
                        end
                        cx = coordinates(i);
                        dx_c = abs(xv - cx); dx_c = min(dx_c, L - dx_c);
                        ix = find(dx_c <= R);
                        if isempty(ix), continue; end
                        ddx = abs(xv(ix) - cx);
                        ddx = ddx - (ddx > L/2) * L;
                        Microstruture(ix(abs(ddx) <= R)) = true;
                    end

                % ---------------------------------------------------------
                case 2
                    nx = ceil(L(1) / resolution);
                    ny = ceil(L(2) / resolution);
                    xv = linspace(0, L(1), nx);
                    yv = linspace(0, L(2), ny);
                    Microstruture = false(nx, ny);

                    for i = 1:n_obj
                        if n_obj > 100 && mod(i, report_step) == 0
                            fprintf('  Progresso: %.1f%%\n', 100*i/n_obj);
                        end
                        cx = coordinates(i,1);
                        cy = coordinates(i,2);
                        dx_c = abs(xv - cx); dx_c = min(dx_c, L(1) - dx_c); ix = find(dx_c <= R);
                        dy_c = abs(yv - cy); dy_c = min(dy_c, L(2) - dy_c); iy = find(dy_c <= R);
                        if isempty(ix) || isempty(iy), continue; end
                        [IX, IY] = ndgrid(ix, iy);
                        ddx = abs(xv(IX) - cx);   ddx = ddx - (ddx > L(1)/2) * L(1);
                        ddy = abs(yv(IY) - cy);   ddy = ddy - (ddy > L(2)/2) * L(2);
                        mask = ddx.^2 + ddy.^2 <= R^2;
                        lin_idx = sub2ind([nx ny], IX(mask), IY(mask));
                        Microstruture(lin_idx) = true;
                    end

                % ---------------------------------------------------------
                case 3
                    nx = ceil(L(1) / resolution);
                    ny = ceil(L(2) / resolution);
                    nz = ceil(L(3) / resolution);
                    xv = linspace(0, L(1), nx);
                    yv = linspace(0, L(2), ny);
                    zv = linspace(0, L(3), nz);
                    Microstruture = false(nx, ny, nz);

                    for i = 1:n_obj
                        if n_obj > 100 && mod(i, report_step) == 0
                            fprintf('  Progresso: %.1f%%\n', 100*i/n_obj);
                        end
                        cx = coordinates(i,1);
                        cy = coordinates(i,2);
                        cz = coordinates(i,3);
                        dx_c = abs(xv - cx); dx_c = min(dx_c, L(1) - dx_c); ix = find(dx_c <= R);
                        dy_c = abs(yv - cy); dy_c = min(dy_c, L(2) - dy_c); iy = find(dy_c <= R);
                        dz_c = abs(zv - cz); dz_c = min(dz_c, L(3) - dz_c); iz = find(dz_c <= R);
                        if isempty(ix) || isempty(iy) || isempty(iz), continue; end
                        [IX, IY, IZ] = ndgrid(ix, iy, iz);
                        ddx = abs(xv(IX) - cx);   ddx = ddx - (ddx > L(1)/2) * L(1);
                        ddy = abs(yv(IY) - cy);   ddy = ddy - (ddy > L(2)/2) * L(2);
                        ddz = abs(zv(IZ) - cz);   ddz = ddz - (ddz > L(3)/2) * L(3);
                        mask = ddx.^2 + ddy.^2 + ddz.^2 <= R^2;
                        lin_idx = sub2ind([nx ny nz], IX(mask), IY(mask), IZ(mask));
                        Microstruture(lin_idx) = true;
                    end

                % ---------------------------------------------------------
                otherwise
                    error('VoxelizeDomain: coordinates must have 1, 2, or 3 columns.');
            end

            fprintf('Voxelizacao concluida.\n');
        end
        function [r_axis, S2_radial, S2_map, phi_vol, std_vol] = ...
                CalcS2Correlation(Microestrutura, resolucao, L)
            %% CalcS2Correlation
            % Compute the isotropic two-point probability function S2(r)
            % for a 3D binary microstructure using the Wiener-Khinchin theorem
            % (FFT-based autocorrelation).
            %
            % This is a direct MATLAB translation of calculate_s2_correlation()
            % from Center2S2func.py, keeping identical numerics.
            %
            % Theory (Torquato, "Random Heterogeneous Materials", 2002):
            %   S2(r) = < I(x) I(x+r) >
            %         = IFFT( |FFT(I)|^2 ) / N_voxels      (Wiener-Khinchin)
            % The 3D map is then radially averaged to give the 1D function S2(r).
            %
            % Syntax:
            %   [r_axis, S2_radial, S2_map, phi_vol, std_vol] = ...
            %       MaterialClass.CalcS2Correlation(Microestrutura, resolucao, L)
            %
            % Inputs:
            %   Microestrutura : logical or numeric 3D array (binary: 0 = matrix,
            %                    1 = solid).  Size must be [nx ny nz].
            %   resolucao      : voxel edge length (isotropic, same units as L).
            %   L              : 1x3 domain size  [Lx  Ly  Lz]  (same units).
            %
            % Outputs:
            %   r_axis    : 1D radial distance vector [resolucao/2 : resolucao : Lmin/2]
            %               in the same units as resolucao / L.
            %   S2_radial : 1D radially averaged S2(r), same length as r_axis.
            %   S2_map    : 3D S2 map before radial averaging (fftshifted, r=0 at centre).
            %   phi_vol   : mean volume fraction  phi = mean(I(:)).
            %   std_vol   : standard deviation of the indicator field.
            %
            % Example:
            %   % single sphere of radius 10 in a 100^3 box, voxel = 1
            %   [nx,ny,nz] = deal(100,100,100);
            %   [X,Y,Z]    = ndgrid(1:nx, 1:ny, 1:nz);
            %   I = sqrt((X-50).^2+(Y-50).^2+(Z-50).^2) <= 10;
            %   [r,S2,~,phi,~] = MaterialClass.CalcS2Correlation(I, 1, [100 100 100]);
            %   figure; plot(r, S2); xline(phi^2,'--'); xlabel('r'); ylabel('S2(r)');
            %
            % See also: MaterialClass.VoxelizeDomain, MaterialClass.GetPSDFromImage

            fprintf('Calculando funcao de correlacao S2(r)...\n');

            % --- global statistics ---
            I_float  = double(Microestrutura);
            phi_vol  = mean(I_float(:));
            std_vol  = std(I_float(:));

            % --- FFT-based autocorrelation (Wiener-Khinchin) ---
            %   S2_map = IFFT( |FFT(I)|^2 ) / N
            if 0
                I_float = I_float - phi_vol;
                I_float = I_float./std_vol;
            end
            F      = fftn(I_float);
            S2_map = real(ifftn(F .* conj(F))) / numel(I_float);

            % shift so that r = 0 is at the array centre
            S2_map = fftshift(S2_map);

            % --- build 3D radial-distance grid (in physical units) ---
            [cx, cy, cz] = size(S2_map);
            center = ([cx cy cz] / 2) - 0.5;   % sub-voxel centre (0-based offset)

            [X_grid, Y_grid, Z_grid] = ndgrid( ...
                (0:cx-1) - center(1), ...
                (0:cy-1) - center(2), ...
                (0:cz-1) - center(3));

            radius_grid = sqrt(X_grid.^2 + Y_grid.^2 + Z_grid.^2) * resolucao;

            % --- radial binning ---
            max_r  = min(L) / 2;
            edges  = 0 : resolucao : (max_r + resolucao);
            r_axis = edges(1:end-1) + resolucao / 2;    % bin centres
            n_bins = numel(r_axis);

            r_flat  = radius_grid(:);
            S2_flat = S2_map(:);

            % assign each voxel to a bin (0-based bin index, then +1 for MATLAB)
            [~, ~, bin_idx] = histcounts(r_flat, edges);

            % keep only voxels that fall inside the range
            valid = bin_idx >= 1 & bin_idx <= n_bins;
            bin_valid = bin_idx(valid);
            S2_valid  = S2_flat(valid);

            % mean S2 per radial bin (accumarray is faster than a loop)
            S2_radial = accumarray(bin_valid, S2_valid, [n_bins 1], @mean, 0);

            fprintf('Calculo de S2(r) concluido.\n');
        end
        function f_hat = hankel0_fwd(f_r, r, k, dr)
            %% hankel0_fwd
            % Forward 2D Hankel (Fourier-Bessel) transform of order 0:
            %   f_hat(k) = 2*pi * int_0^inf f(r) * J0(k*r) * r  dr
            %
            % Used internally by MonoDisperseSphere for the 2D PY OZ solver.
            % Block-wise direct quadrature keeps memory usage bounded.
            %
            % Inputs:
            %   f_r : column vector of f(r) values on grid r
            %   r   : real-space grid   (column vector)
            %   k   : wavenumber grid   (column vector)
            %   dr  : real-space step size
            %
            % Output:
            %   f_hat : column vector of transform values on grid k
            Nk    = length(k);
            w     = 2*pi * r * dr;          % quadrature weights (Nr x 1)
            fw    = f_r .* w;               % pre-weighted integrand
            f_hat = zeros(Nk, 1);
            block = 256;
            for i0 = 1:block:Nk
                i1 = min(i0 + block - 1, Nk);
                J0 = besselj(0, k(i0:i1) * r');   % (blk x Nr)
                f_hat(i0:i1) = J0 * fw;
            end
        end
        function f_r = hankel0_inv(f_k, k, r, dk)
            %% hankel0_inv
            % Inverse 2D Hankel transform of order 0:
            %   f(r) = (1/(2*pi)) * int_0^inf f_hat(k) * J0(k*r) * k  dk
            %
            % Used internally by MonoDisperseSphere for the 2D PY OZ solver.
            %
            % Inputs:
            %   f_k : column vector of f_hat(k) values on grid k
            %   k   : wavenumber grid   (column vector)
            %   r   : real-space grid   (column vector)
            %   dk  : wavenumber step size
            %
            % Output:
            %   f_r : column vector of f(r) values on grid r
            Nr  = length(r);
            w   = (1/(2*pi)) * k * dk;     % quadrature weights (Nk x 1)
            fw  = f_k .* w;                % pre-weighted integrand
            f_r = zeros(Nr, 1);
            block = 256;
            for i0 = 1:block:Nr
                i1 = min(i0 + block - 1, Nr);
                J0 = besselj(0, r(i0:i1) * k');   % (blk x Nk)
                f_r(i0:i1) = J0 * fw;
            end
        end
        function [g_r, S_k, r_grid, k_grid] = HardBodyPY(eta, D, d)
            %% HardBodyPY
            % Pair correlation function g(r) and static structure factor S(k)
            % for monodisperse hard bodies in d dimensions via the Percus-Yevick
            % (PY) Ornstein-Zernike (OZ) integral equation.
            %
            %   d = 1 : Hard rods     — 1D OZ+PY Picard iteration, cosine transforms
            %   d = 2 : Hard disks    — 2D OZ+PY Picard iteration, Hankel-0 transforms
            %   d = 3 : Hard spheres  — analytic 3D PY c(k), sine transform for g(r)
            %
            % Syntax:
            %   [g_r, S_k, r_grid, k_grid] = MaterialClass.HardBodyPY(eta, D, d)
            %
            % Inputs:
            %   eta    – packing fraction (length / area / volume fraction)
            %   D      – object diameter (or rod length for d=1)
            %   d      – spatial dimension (1, 2, or 3)
            %
            % Outputs:
            %   g_r    – pair correlation function values on r_grid
            %   S_k    – static structure factor values on k_grid
            %   r_grid – physical r values (same units as D)
            %   k_grid – wavenumber values; physical (1/D) for d=1,2;
            %            k_norm = k_phys*(D/2) for d=3 (matches MonoDisperseSphere)

            switch d
                case 1
                    % ---- 1D hard rods: Picard iteration with cosine transforms ----
                    Nr    = 4096;
                    r_max = 20.0 * D;
                    dr    = r_max / Nr;
                    r_grid = ((1:Nr) - 0.5)' * dr;

                    dk = pi / r_max;
                    Nk = floor((pi / dr) / dk);
                    k_grid = ((1:Nk) - 0.5)' * dk;

                    rho_1D  = eta / D;
                    inside  = r_grid < D;
                    outside = ~inside;
                    alpha   = 0.4;
                    tol     = 1e-12;

                    cos_fwd = @(f) MaterialClass.cosine_fwd_1D(f, r_grid, k_grid, dr);
                    cos_inv = @(f) MaterialClass.cosine_inv_1D(f, k_grid, r_grid, dk);

                    gamma_r = zeros(Nr, 1);
                    c_r = zeros(Nr, 1);
                    c_r(inside) = -1.0;

                    converged = false;
                    for iter = 1:5000
                        c_hat = cos_fwd(c_r);
                        h_hat = c_hat ./ (1.0 - rho_1D * c_hat);
                        h_r   = cos_inv(h_hat);
                        gamma_r_new = h_r - c_r;
                        err = max(abs(gamma_r_new - gamma_r));
                        gamma_r = (1 - alpha)*gamma_r + alpha*gamma_r_new;
                        c_r(inside)  = -1.0 - gamma_r(inside);
                        c_r(outside) = 0.0;
                        if err < tol
                            fprintf('    1D PY converged at iter %d, err = %.3e\n', iter, err);
                            converged = true;  break;
                        end
                    end
                    if ~converged
                        fprintf('    WARNING: 1D PY did not converge (err=%.3e).\n', err);
                    end

                    c_hat = cos_fwd(c_r);
                    h_hat = c_hat ./ (1.0 - rho_1D * c_hat);
                    h_r   = cos_inv(h_hat);
                    g_r   = 1.0 + h_r;
                    g_r(inside) = 0.0;
                    S_k   = 1.0 + rho_1D * h_hat;

                case 2
                    % ---- 2D hard disks: Picard iteration with Hankel-0 transforms ----
                    Nr    = 4096;
                    r_max = 20.0 * D;
                    dr    = r_max / Nr;
                    r_grid = ((1:Nr) - 0.5)' * dr;

                    dk = pi / r_max;
                    Nk = floor((pi / dr) / dk);
                    k_grid = ((1:Nk) - 0.5)' * dk;

                    rho_2D  = eta / (pi*(D/2)^2);
                    inside  = r_grid < D;
                    outside = ~inside;
                    alpha   = 0.4;
                    tol     = 1e-12;

                    h_fwd = @(f) MaterialClass.hankel0_fwd(f, r_grid, k_grid, dr);
                    h_inv = @(f) MaterialClass.hankel0_inv(f, k_grid, r_grid, dk);

                    gamma_r = zeros(Nr, 1);
                    c_r = zeros(Nr, 1);
                    c_r(inside) = -1.0;

                    converged = false;
                    for iter = 1:5000
                        c_k = h_fwd(c_r);
                        h_k = c_k ./ (1.0 - rho_2D * c_k);
                        h_r = h_inv(h_k);
                        gamma_r_new = h_r - c_r;
                        err = max(abs(gamma_r_new - gamma_r));
                        gamma_r = (1 - alpha)*gamma_r + alpha*gamma_r_new;
                        c_r(inside)  = -1.0 - gamma_r(inside);
                        c_r(outside) = 0.0;
                        if err < tol
                            fprintf('    2D PY converged at iter %d, err = %.3e\n', iter, err);
                            converged = true;  break;
                        end
                    end
                    if ~converged
                        fprintf('    WARNING: 2D PY did not converge (err=%.3e).\n', err);
                    end

                    c_k = h_fwd(c_r);
                    h_k = c_k ./ (1.0 - rho_2D * c_k);
                    h_r = h_inv(h_k);
                    g_r = 1.0 + h_r;
                    g_r(inside) = 0.0;
                    S_k = 1.0 + rho_2D * h_k;

                case 3
                    % ---- 3D hard spheres: analytic PY c(k) + sine transform for g(r) ----
                    % k_grid stores k_norm = k_phys*(D/2), range 0 to ~3
                    % r_grid stores physical r, range 0 to 5*D
                    k_grid = (linspace(0, 6/D, 4094))';
                    r_grid = (linspace(0, 5*D,  4096))';

                    k_norm = k_grid * (D/2);   % dimensionless k_norm = k_phys*(D/2)
                    r_norm = r_grid * (2/D);    % dimensionless r_norm = r_phys/(D/2)

                    rhoS = 3*eta / (4*pi);
                    l1   = (1 + 2*eta)^2 / (1 - eta)^4;
                    l2   = -(1 + eta/2)^2 / (1 - eta)^4;

                    % Analytic PY direct correlation c_hat(k_norm)
                    c_hat = -4*pi ./ (k_norm.^3) .* ( ...
                        l1*(sin(2*k_norm) - 2*k_norm.*cos(2*k_norm)) + ...
                        3*eta*l2./k_norm .* (4*k_norm.*sin(2*k_norm) + (2 - 4*k_norm.^2).*cos(2*k_norm) - 2) + ...
                        eta*l1./(2*k_norm.^3) .* ((6*k_norm.^2 - 3 - 2*k_norm.^4).*cos(2*k_norm) + ...
                        (4*k_norm.^3 - 6*k_norm).*sin(2*k_norm) + 3) );
                    c_hat(k_norm < 1e-8) = -8*pi/3 * ((4 + eta)*l1 + 18*eta*l2);

                    h3D = c_hat ./ (1 - rhoS * c_hat);   % OZ relation
                    S_k = 1 + rhoS * h3D;                 % structure factor

                    % g(r) via inverse 3D sine transform:
                    %   h(r_norm) = 1/(2*pi^2*r_norm) * int h_hat(k_norm)*k_norm*sin(k*r) dk
                    dk_norm = k_norm(2) - k_norm(1);
                    h_r = zeros(numel(r_norm), 1);
                    % r_norm = 0: use l'Hopital limit sin(k*r)/r → k
                    h_r(1) = trapz(k_norm, h3D .* k_norm.^2) / (2*pi^2);
                    for ii = 2:numel(r_norm)
                        h_r(ii) = trapz(k_norm, h3D .* k_norm .* sin(k_norm * r_norm(ii))) ...
                                  / (2*pi^2 * r_norm(ii));
                    end
                    g_r = max(0, 1 + h_r);
                    g_r(r_grid < D) = 0;   % enforce hard-core exclusion

                otherwise
                    error('HardBodyPY: d must be 1, 2, or 3.');
            end
        end
        function f_hat = cosine_fwd_1D(f_r, r, k, dr)
            %% cosine_fwd_1D
            % Forward 1D cosine (Fourier) transform for even functions:
            %   f_hat(k) = 2 * int_0^inf f(r) cos(kr) dr
            %
            % Used internally by HardBodyPY for the 1D PY OZ solver.
            %
            % Inputs:
            %   f_r : column vector of f(r) values on grid r
            %   r   : real-space grid   (column vector)
            %   k   : wavenumber grid   (column vector)
            %   dr  : real-space step size
            %
            % Output:
            %   f_hat : column vector of transform values on grid k
            Nk    = length(k);
            w     = 2 * dr;
            fw    = f_r * w;
            f_hat = zeros(Nk, 1);
            block = 256;
            for i0 = 1:block:Nk
                i1 = min(i0 + block - 1, Nk);
                C = cos(k(i0:i1) * r');   % (blk x Nr)
                f_hat(i0:i1) = C * fw;
            end
        end
        function f_r = cosine_inv_1D(f_hat, k, r, dk)
            %% cosine_inv_1D
            % Inverse 1D cosine (Fourier) transform for even functions:
            %   f(r) = (1/pi) * int_0^inf f_hat(k) cos(kr) dk
            %
            % Used internally by HardBodyPY and MonoDisperseSphere (d=1).
            %
            % Inputs:
            %   f_hat : column vector of f_hat(k) values on grid k
            %   k     : wavenumber grid   (column vector)
            %   r     : real-space grid   (column vector)
            %   dk    : wavenumber step size
            %
            % Output:
            %   f_r : column vector of f(r) values on grid r
            Nr  = length(r);
            w   = dk / pi;
            fw  = f_hat * w;
            f_r = zeros(Nr, 1);
            block = 256;
            for i0 = 1:block:Nr
                i1 = min(i0 + block - 1, Nr);
                C = cos(r(i0:i1) * k');   % (blk x Nk)
                f_r(i0:i1) = C * fw;
            end
        end
        function mat = prepareSigma(mat,d)
            %% prepareSigma
            % Prepare scattering cross-sections and derived quantities
            %
            % Syntax:
            %   mat = MaterialClass.prepareSigma( mat, d );
            %
            % Inputs:
            %   mat : MaterialClass object
            %   d   : dimension of the problem

            if isempty(mat.sigma)
                mat.CalcSigma;
            end

            if mat.acoustics
                [mat.Sigma,mat.Sigmapr,mat.invcdf] = MaterialClass.prepareSigmaOne(mat.sigma{1},d);

                if mat.Sigma == 0
                    % Homogeneous acoustic medium: propagation is ballistic.
                    mat.meanFreeTime = Inf;
                    mat.meanFreePath = Inf;
                    mat.transportMeanFreeTime = Inf;
                    mat.transportMeanFreePath = Inf;

                    % Diffusivity and scattering anisotropy are undefined
                    % when no normalized scattering phase function exists.
                    mat.Diffusivity = NaN;
                    mat.g = NaN;
                else
                    % Diffusion coefficient m²/s (Eq. (5.12), Ryzhik et al, 1996)
                    mat.Diffusivity = mat.v^2/(double(d)*(mat.Sigma-mat.Sigmapr));

                    mat.meanFreeTime = 1/mat.Sigma;
                    mat.meanFreePath = mat.v * mat.meanFreeTime;

                    mat.transportMeanFreePath = double(d) * mat.Diffusivity / mat.v;
                    mat.transportMeanFreeTime = mat.transportMeanFreePath/mat.v;

                    % anisotropy coefficient (characterizes scattering directionality)
                    mat.g = 1 - mat.meanFreePath/mat.transportMeanFreePath;
                end

                % the two lines below are just for homogenization of the propagation
                % code between acoustics and elastics
                mat.vp = mat.v;
                mat.vs = 0;
            else
                mat.Sigma = zeros(2);
                mat.Sigmapr = zeros(2);
                mat.invcdf = cell(2);

                [mat.Sigma(1,1),mat.Sigmapr(1,1),mat.invcdf{1,1}] = MaterialClass.prepareSigmaOne(mat.sigma{1,1},d);
                [mat.Sigma(1,2),mat.Sigmapr(1,2),mat.invcdf{1,2}] = MaterialClass.prepareSigmaOne(mat.sigma{1,2},d);
                [mat.Sigma(2,1),mat.Sigmapr(2,1),mat.invcdf{2,1}] = MaterialClass.prepareSigmaOne(mat.sigma{2,1},d);
                [mat.Sigma(2,2),mat.Sigmapr(2,2),mat.invcdf{2,2}] = MaterialClass.prepareSigmaOne(mat.sigma{2,2},d);
                SigmaP = mat.Sigma(1,1) + mat.Sigma(1,2);
                SigmaS = mat.Sigma(2,1) + mat.Sigma(2,2);
                mat.meanFreeTime = 1./[SigmaP SigmaS];
                mat.meanFreePath = [mat.vp mat.vs].' .* mat.meanFreeTime;

                % Safe same-mode probabilities for rows having no scattering
                mat.P2P = 1;
                mat.S2S = 1;
                if SigmaP > 0
                    mat.P2P = mat.Sigma(1,1)/SigmaP;
                end
                if SigmaS > 0
                    mat.S2S = mat.Sigma(2,2)/SigmaS;
                end

                if all([SigmaP SigmaS] == 0)
                    % Homogeneous elastic medium: neither mode scatters
                    mat.transportMeanFreePath = [Inf; Inf];
                    mat.transportMeanFreeTime = [Inf; Inf];
                    mat.Diffusivity = NaN;
                elseif any([SigmaP SigmaS] == 0)
                    % One ballistic mode prevents use of the coupled
                    % two-mode diffusion approximation
                    mat.transportMeanFreePath = [NaN; NaN];
                    mat.transportMeanFreeTime = [NaN; NaN];
                    mat.Diffusivity = NaN;
                else
                    K = mat.vp/mat.vs;
                    % Transport (diffusion) mean free paths of P & S waves
                    tmfp_P = (mat.vp*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2)) + mat.vs*mat.Sigmapr(1,2) )...
                                      /( (mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1))*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2))- mat.Sigmapr(1,2)*mat.Sigmapr(2,1) );
                    tmfp_S = (mat.vs*(mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1)) + mat.vp*mat.Sigmapr(2,1) )...
                                      /( (mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1))*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2))- mat.Sigmapr(1,2)*mat.Sigmapr(2,1) );
                    mat.transportMeanFreePath = [tmfp_P tmfp_S]';
                    mat.transportMeanFreeTime = mat.transportMeanFreePath ./ ...
                        [mat.vp mat.vs].';

                    % partial diffusion coefficients of P & S waves
                    Dp = mat.vp*tmfp_P/d;
                    Ds = mat.vs*tmfp_S/d;
                    % Diffusion coefficient m²/s (Eqs. (5.42) & (5.46), Ryzhik et al, 1996)
                    mat.Diffusivity = double((Dp+2*K^3*Ds)/(1+2*K^3));
                end
            end
        end
        function [Sigma,Sigma_prime,invcdf] = prepareSigmaOne(sigma,d)
            if d==2
                Sigma = 2*integral(sigma,0,pi);
                Sigma_prime = 2*integral(@(th)sigma(th).*cos(th),0,pi);
            elseif d==3
                Sigma = 2*pi*integral(@(th)sigma(th).*sin(th),0,pi);
                Sigma_prime = 2*pi*integral(@(th)sigma(th).*sin(th).*cos(th),0,pi);
            else
                error('MaterialClass:prepareSigmaOne:InvalidDimension', ...
                    'The scattering preparation requires d = 2 or d = 3.');
            end

            if ~isfinite(Sigma) || Sigma < 0
                error('MaterialClass:prepareSigmaOne:InvalidSigma', ...
                    'The integrated scattering cross-section must be finite and nonnegative.');
            end

            if Sigma == 0
                % There is no angular probability distribution to invert.
                % This placeholder is never sampled because the associated
                % mean free time is infinite.
                Sigma_prime = 0;
                invcdf = @(probability) zeros(size(probability));
                return
            end

            if d==2
                sigmaNorm = @(th) (2/Sigma)*sigma(th);
            else
                sigmaNorm = @(th) (2*pi/Sigma)*sin(th).*sigma(th);
            end

            Nth = 1e6;
            xth = linspace(0,pi,Nth);
            pdf = sigmaNorm(xth);
            if any(isnan(pdf))
                warning('There is a NaN inside the probability density function')
            end
            % Build the scattering-angle CDF.
            cdf = cumsum(pdf,"omitnan")*mean(diff(xth));

            % Normalize the scattering-angle CDF.
            cdfTotal = cdf(end);
            if ~isfinite(cdfTotal) || cdfTotal <= 0
                error('MaterialClass:prepareSigmaOne:InvalidCDF', ...
                    'Could not build a valid scattering-angle CDF.');
            end
            cdf = cdf ./ cdfTotal;
            cdf(1) = 0;
            cdf(end) = 1;

            % Remove repeated CDF values before constructing its inverse.
            ind = find(diff(cdf)>0);
            ind = unique([1 ind ind+1 Nth]);
            [cdfUnique, idx] = unique(cdf(ind),'stable');
            xthUnique = xth(ind);
            invcdf = griddedInterpolant(cdfUnique,xthUnique(idx));
        end
        function centers = CreateSphereComposite(L,D,phi)
            %cria um composito de matriz m (prop1) e inclusao i (prop2)
            %usando o metodo Random Sequential 

            %L is the size can be 2D or 3D
            % D is the diameter
            % phi is the volume fraction

            dim = numel(L);

            switch dim
                case 1
                    [centers, ~] = MaterialClass.rod_packing_ls(L,D,phi);
                case 2
                    [centers, ~] = MaterialClass.disk_packing_ls(L,D,phi);
                case 3
                    [centers, ~] = MaterialClass.sphere_packing_ls(L,D,phi);
            end
        end
        [centers, nobj] = rod_packing_ls(L,D,phi);
        [centers, nobj] = disk_packing_ls(L,D,phi);
        [centers, nobj] = sphere_packing_ls(L,D,phi);
        function psd_summary_figure(dim_label, phi_vol, Lc, r_plot, S2_plot, R_plot, k_vec, psd_vals, phi)
            figure('Name', sprintf('GetPSDFromImage — %s Microstructure', dim_label), ...
                'Color', 'w', 'Position', [120 80 1200 400]);
            tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
            title(tl, sprintf('%s Image PSDF  (\\phi = %.4f,  L_c = %.4g)', ...
                dim_label, phi_vol, Lc), 'FontSize', 13, 'FontWeight', 'bold');

            nexttile;
            plot(r_plot, S2_plot, 'm-', 'LineWidth', 1.8);
            hold on;
            yline(phi^2, 'k--', 'LineWidth', 0.8, 'Label', '\phi^2', ...
                'LabelVerticalAlignment', 'bottom');
            yline(phi, 'k:', 'LineWidth', 0.8, 'Label', '\phi');
            xlabel('r / L_c');  ylabel('S_2(r)');
            title('Two-Point Correlation S_2(r)');
            xlim([0 max(r_plot)]);  grid on;  box on;

            nexttile;
            plot(r_plot, R_plot, 'b-', 'LineWidth', 1.8);
            hold on;
            yline(0, 'k--', 'LineWidth', 0.8);
            xlabel('r / L_c');  ylabel('R(r)');
            title('Normalised Autocorrelation R(r)');
            xlim([0 max(r_plot)]);  grid on;  box on;

            nexttile;
            plot(k_vec, psd_vals, 'r-', 'LineWidth', 1.8);
            xlabel('k \cdot L_c');  ylabel('\Phi(k)');
            title('Power Spectral Density \Phi(k)');
            xlim([0 min(6, k_vec(end))]);  grid on;  box on;
        end
        function plot_map(centers, DD, L)
            %% plot_map  Visualise a packing produced by the *_packing_ls routines.
            %
            % Syntax:
            %   MaterialClass.plot_map(centers, DD, L)
            %
            % Inputs:
            %   centers  – n-by-d array of object centres (d = 1, 2, or 3)
            %   DD       – object diameter (or length for 1D)
            %   L        – box size: scalar (1D), 1x2 (2D), or 1x3 (3D)

            dim = size(centers, 2);
            n = size(centers,1);
            switch dim
                case 1
                    figure('Name', '1D Rod Packing', 'Color', 'w');
                    hold on;
                    ylim([0 2]);
                    xlim([0 L]);
                    title(sprintf('1D Hard Rod Packing — N=%d, D=%.4g, L=%.4g, \\phi=%.4f', ...
                        n, DD, L, n*DD/L));
                    for i = 1:n
                        x = centers(i);
                        rectangle('Position', [x-DD/2, 0.8, DD, 0.4], ...
                            'FaceColor', [0.3 0.6 1.0], 'EdgeColor', [0.1 0.3 0.7], 'LineWidth', 1.5);
                        if x - DD/2 < 0
                            rectangle('Position', [x-DD/2+L, 0.8, DD-(x-DD/2+L-L), 0.4], ...
                                'FaceColor', [0.3 0.6 1.0], 'EdgeColor', [0.1 0.3 0.7]);
                        end
                        if x + DD/2 > L
                            rectangle('Position', [0, 0.8, x+DD/2-L, 0.4], ...
                                'FaceColor', [0.3 0.6 1.0], 'EdgeColor', [0.1 0.3 0.7]);
                        end
                    end
                    line([0 L], [0.7 0.7], 'Color', 'k', 'LineWidth', 2);
                    xlabel('Position');
                    set(gca, 'YTick', []);
                    hold off;

                case 2
                    figure('Name', 'Disk Packing', 'Color', 'w');
                    hold on; axis equal;
                    xlim([0 L(1)]); ylim([0 L(2)]);
                    xlabel('x'); ylabel('y');
                    phi = n * pi * (DD/2)^2 / (L(1)*L(2));
                    title(sprintf('Disk Packing — N=%d, D=%.4g, Box=%.2g x %.2g, \\phi=%.4f', ...
                        n, DD, L(1), L(2), phi));
                    theta = linspace(0, 2*pi, 60);
                    r  = DD / 2;
                    cx = r * cos(theta);
                    cy = r * sin(theta);
                    for i = 1:n
                        fill(centers(i,1)+cx, centers(i,2)+cy, [0.3 0.6 1.0], ...
                            'EdgeColor', [0.1 0.3 0.7], 'LineWidth', 0.4, 'FaceAlpha', 0.7);
                    end
                    rectangle('Position', [0, 0, L(1), L(2)], 'EdgeColor', 'k', 'LineWidth', 1.5);
                    hold off;

                case 3
                    figure('Name', 'Sphere Packing', 'Color', 'w');
                    hold on; axis equal; grid on;
                    xlim([0 L(1)]); ylim([0 L(2)]); zlim([0 L(3)]);
                    xlabel('x'); ylabel('y'); zlabel('z');
                    phi = n * (4/3)*pi*(DD/2)^3 / (L(1)*L(2)*L(3));
                    title(sprintf('Sphere Packing — N=%d, D=%.4g, \\phi=%.4f', n, DD, phi));
                    marker_size = (DD / max(L))*100;
                    scatter3(centers(1:n,1), centers(1:n,2), centers(1:n,3), ...
                        marker_size, [0.3 0.6 1.0], 'filled', ...
                        'MarkerEdgeColor', [0.1 0.3 0.7], 'MarkerFaceAlpha', 0.7);
                    corners = [0 0 0; L(1) 0 0; L(1) L(2) 0; 0 L(2) 0; ...
                               0 0 L(3); L(1) 0 L(3); L(1) L(2) L(3); 0 L(2) L(3)];
                    edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
                    for e = 1:size(edges, 1)
                        line(corners(edges(e,:),1)', corners(edges(e,:),2)', corners(edges(e,:),3)', ...
                            'Color', 'k', 'LineWidth', 1.5);
                    end
                    view(3);
                    hold off;

                otherwise
                    error('plot_map: centers must have 1, 2, or 3 columns.');
            end
        end
        [med_lam, std_lam, med_mu, std_mu, med_rho, std_rho] = convert_kmr_to_lmr(med_k, std_k, med_mu_in, std_mu_in, med_rho_in, std_rho_in,corr_K_mu)

        function h = PlotVoxel(M, resolution, slice_index)
            %% PlotVoxel  Visualise a voxelised/pixelised microstructure.
            %
            % For a 2D binary array the full image is shown.
            % For a 3D binary volume an XY cross-section is shown.
            %
            % Syntax:
            %   MaterialClass.PlotVoxel(M, resolution)
            %   MaterialClass.PlotVoxel(M, resolution, slice_index)
            %
            % Inputs:
            %   M           : logical/numeric array — nx×ny (2D) or nx×ny×nz (3D)
            %   resolution  : voxel/pixel edge length (isotropic, same unit as domain)
            %   slice_index : z-index of the XY slice to display (3D only).
            %                 Defaults to the middle slice.

            if nargin < 2 || isempty(resolution), resolution = 1; end

            nd = ndims(M);
            if nd == 2 || (nd == 3 && size(M, 3) == 1)
                % ---- 2D image ----
                [nx, ny] = size(M);
                x_ax = (0:ny-1) * resolution;
                y_ax = (0:nx-1) * resolution;
                h = figure('Name', '2D Microstructure', 'Color', 'w');
                imagesc(x_ax, y_ax, double(M));
                colormap(gray(2));
                clim([0 1]);
                axis image;
                xlabel('x'); ylabel('y');
                phi = mean(M(:));
                title(sprintf('2D microstructure  (\\phi = %.4f)', phi));
                colorbar('Ticks', [0.25 0.75], 'TickLabels', {'matrix','inclusion'});
            else
                % ---- 3D volume — one XY slice ----
                [nx, ny, nz] = size(M);
                if nargin < 3 || isempty(slice_index)
                    slice_index = round(nz / 2);
                end
                slice_index = max(1, min(nz, slice_index));

                x_ax = (0:ny-1) * resolution;
                y_ax = (0:nx-1) * resolution;
                h = figure('Name', sprintf('3D Microstructure — XY slice z=%d/%d', slice_index, nz), ...
                    'Color', 'w');
                imagesc(x_ax, y_ax, double(M(:,:,slice_index)));
                colormap(gray(2));
                clim([0 1]);
                axis image;
                xlabel('x'); ylabel('y');
                phi = mean(M(:));
                title(sprintf('3D microstructure — XY slice z=%d/%d  (\\phi_{vol} = %.4f)', ...
                    slice_index, nz, phi));
                colorbar('Ticks', [0.25 0.75], 'TickLabels', {'matrix','inclusion'});
            end
            set(gca, 'FontSize', 13);
            box on;
        end
    end
end
