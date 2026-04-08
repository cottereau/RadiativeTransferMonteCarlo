%% MaterialClass
% Class to deal with Differential Scattering Cross-Sections (DSCS)
%
%
% MaterialClass ()
%
% See also
%
%
classdef MaterialClass < handle
    properties
        % Required Properties
        d               int8    = 3; % dimension
        %mat             struct  = struct.empty % PSDF parameter
        type            char    = 'isotropic'
        acoustics

        v
        vp
        vs
        rho
        Frequency
        coefficients_of_variation
        correlation_coefficients

        SpectralLaw     char    = ''; % PSDF name
        SpectralParam   struct  = struct.empty % PSDF parameter

        CorrelationLength       = []; % correlation length


        sigma           cell    = cell.empty; % Differential Scattering Cross-Sections

        Sigma
        Sigmapr
        invcdf
        Diffusivity             double = [];
        meanFreeTime            double = [];
        meanFreePath            double = [];
        transportMeanFreeTime   double = [];
        transportMeanFreePath   double = [];
        g                       double = [];
        P2P
        S2S

        Phi                     = []; % power spectral density / function_handle
        k               double  = []; % wavenumber vector
        R                       = []; % Correlation / function_handle
        r               double  = []; % r vector
        timeSteps               =  0; % time Steps : 0=small 1=large

    end
    properties (SetAccess = private, Hidden = true)
        Type_def = {'isotropic'}; %the anisotropic should be implemented
        SpectralLaw_def = {'exp','power_law','gaussian','triangular','low_pass','VonKarman','monodispersesphere','image','Imported'};
    end
    methods
        function obj = MaterialClass(geometry,freq,acoustics, ...
                v,coefficients_of_variation,correlation_coefficients, ...
                acf,lc)
            %% MaterialClass
            % MaterialClass contructor
            %
            % Syntax:
            %   newobj = MaterialClass (  );
            %
            % Inputs:
            %     mat : material structure
            %     d   : dimension of the problem

            if nargin ~=0
                if geometry.dimension~=3
                    warning(['Total scattering cross-sections for 2D case are to be ' ...
                        'implemented in future releases of the code'])
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
        %% GET AND SET METHODS
        function set.type(obj,newvalue)
            obj.type = validatestring(newvalue,obj.Type_def); %#ok<MCSUP>
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

            % get the power spectral density function
            if isempty(obj.Phi)
                obj.getPSDF;
            end

            switch obj.type
                case 'isotropic'
                    obj.calcSigmaIsotropic;
                case 'anisotropic'
                    obj.calcSigmaAnisotropic;
                otherwise
                    error('DSCS type not supported')
            end
        end
        function calcSigmaIsotropic(obj)
            %% calcSigmaIsotropic
            % Compute the normalized differential scattering cross-sections
            % based on statistics of the fluctuating material parameters
            % of the wave equation for isotropic PSDF
            %
            % Syntax:
            %   newobj = calcSigmaIsotropic( );
            %
            % Inputs:
            %
            % Output:
            % The 3D formulas are based on:
            %   L. Ryzhik, G. Papanicolaou, J. B. Keller. Transport equations for elastic
            %   and other waves in random media. Wave Motion 24, pp. 327-370, 1996.
            %   doi: 10.1016/S0165-2125(96)00021-2
            %
            % The 2D formulas are based on:
            %   H. Sato, M. C. Fehler, T. Maeda. Seismic Wave Propagation and Scattering
            %   in the Heterogeneous Earth, 2nd edition. Springer, 2012.
            %   (Chapter 4, Born approximation for wave scattering in random media)
            %   See also: R. S. Wu and K. Aki. Scattering characteristics of elastic waves
            %   by an elastic heterogeneity. Geophysics 50(4), pp. 582-595, 1985.
            %   doi: 10.1190/1.1441934

            switch obj.acoustics
                case 1
                    switch obj.d
                        case 1
                            error('\sigma for 1D not implemented')
                        case 2
                            % 2D ACOUSTICS
                            % In 2D the differential scattering cross section has units of
                            % length (not area), and the solid angle is replaced by dtheta.
                            % Leading prefactor is (pi/4)*zeta^2 instead of (pi/2)*zeta^3.
                            % Reference: Sato, Fehler & Maeda (2012), Ch. 4
                            %
                            %   sigma(th) = (pi/4) * zeta^2 * omega
                            %               * (cos(th)^2 * delta_rr^2
                            %                  + 2*cos(th)*rho_kr*delta_kk*delta_rr
                            %                  + delta_kk^2)
                            %               * Phi(zeta * sqrt(2*(1-cos(th))))
                            %
                            % This is the 2D analogue of Ryzhik et al. (1996) Eq. (1.3).
                            % Note: Phi is the 2D isotropic PSDF evaluated at the
                            % scalar wavenumber |k' - k| = k*sqrt(2*(1-cos(th))).

                            zeta = (2*pi*obj.Frequency)/obj.v * obj.CorrelationLength;
                            delta_kk = obj.coefficients_of_variation(1); % CV of compressibility
                            delta_rr = obj.coefficients_of_variation(2); % CV of density
                            rho_kr   = obj.correlation_coefficients;     % corr. of kappa and rho

                            obj.sigma = {@(th) (pi/4)*zeta^2 * ...
                                (cos(th).^2*delta_rr^2 + ...
                                2*cos(th)*rho_kr*delta_kk*delta_rr + delta_kk^2) ...
                                .*obj.Phi(zeta.*sqrt(2*(1-cos(th)))) * (2*pi*obj.Frequency)};

                        case 3
                            % acoustics
                            zeta = (2*pi*obj.Frequency)/obj.v*obj.CorrelationLength;
                            delta_kk = obj.coefficients_of_variation(1); % coefficient of variation of compressibility
                            delta_rr = obj.coefficients_of_variation(2); % coefficient of variation of density
                            rho_kr = obj.correlation_coefficients; % correlation of compressibility and density

                            % [Ryzhik et al, 1996; Eq. (1.3)]
                            obj.sigma = {@(th) (pi/2)*zeta^3*(cos(th).^2*delta_rr^2 + ...
                                2*cos(th)*rho_kr*delta_kk*delta_rr + delta_kk^2) ...
                                .*obj.Phi(zeta.*sqrt(2*(1-cos(th))))*(2*pi*obj.Frequency)};
                    end
                case 0
                    switch obj.d
                        case 1
                            error('\sigma for 1D not implemented')
                        case 2
                            % 2D ELASTICS
                            % In 2D the prefactor changes from (pi/2)*zeta^3 -> (pi/4)*zeta^2,
                            % and similarly for S-wave terms.
                            % Reference: Sato, Fehler & Maeda (2012), Ch. 4;
                            %            Wu & Aki (1985), Geophysics 50(4).
                            %
                            % All angular radiation-pattern numerators are identical to the
                            % 3D case (they come from the polarisation/geometry of the
                            % incident and scattered modes); only the dimensional prefactor
                            % and the PSDF normalisation change.

                            K      = obj.vp / obj.vs;
                            zetaP  = (2*pi*obj.Frequency) / obj.vp * obj.CorrelationLength;
                            zetaS  = K * zetaP;

                            delta_ll = obj.coefficients_of_variation(1); % CV of lambda
                            delta_mm = obj.coefficients_of_variation(2); % CV of mu
                            delta_rr = obj.coefficients_of_variation(3); % CV of density
                            rho_lm   = obj.correlation_coefficients(1);  % corr. lambda-mu
                            rho_lr   = obj.correlation_coefficients(2);  % corr. lambda-density
                            rho_mr   = obj.correlation_coefficients(3);  % corr. mu-density

                            % --- P -> P ---
                            % 2D analogue of [Ryzhik et al. 1996, Eq. (1.3)] and
                            % [Sato et al. 2012, Ch. 4 / Wu & Aki 1985]
                            sigmaPP = @(th) (pi/4)*zetaP^2 * ...
                                ( (1-2/K^2)^2*delta_ll^2 ...
                                + 4*(1/K^2 - 2/K^4)*rho_lm*delta_ll*delta_mm.*cos(th).^2 ...
                                + (4/K^4)*delta_mm^2.*cos(th).^4 ...
                                + delta_rr^2.*cos(th).^2 ...
                                + 2*(1-2/K^2)*rho_lr*delta_ll*delta_rr.*cos(th) ...
                                + (4/K^2)*rho_mr*delta_mm*delta_rr.*cos(th).^3 ) ...
                                .*obj.Phi(zetaP.*sqrt(2*(1-cos(th)))) * (2*pi*obj.Frequency);

                            % --- P -> S ---
                            % 2D analogue of [Sato et al. 2012, Ch. 4; Ryzhik et al. Eqs.(4.56),(1.20),(1.22)]
                            sigmaPS = @(th) (pi/4)*K*zetaP^2 * ...
                                ( K^2*delta_rr^2 ...
                                + 4*delta_mm^2.*cos(th).^2 ...
                                + 4*K*rho_mr*delta_mm*delta_rr.*cos(th) ) ...
                                .*(1-cos(th).^2) ...
                                .*obj.Phi(zetaP.*sqrt(1+K^2-2*K*cos(th))) * (2*pi*obj.Frequency);

                            % --- S -> P ---
                            % 2D analogue of [Sato et al. 2012, Ch. 4; Ryzhik et al. Eqs.(4.56),(1.20),(1.21)]
                            sigmaSP = @(th) (pi/8/K^3)*zetaS^2 * ...
                                ( delta_rr^2 ...
                                + (4/K^2)*delta_mm^2.*cos(th).^2 ...
                                + (4/K)*rho_mr*delta_mm*delta_rr.*cos(th) ) ...
                                .*(1-cos(th).^2) ...
                                .*obj.Phi(zetaS.*sqrt(1+1/K^2-2/K*cos(th))) * (2*pi*obj.Frequency);

                            % --- S -> S ---
                            % 2D analogue of [Ryzhik et al., Eq. (4.54); Sato et al. 2012, Ch. 4]
                            % In 2D, SH and SV are decoupled; the combined SS term below
                            % corresponds to the in-plane (SV-SV) scattering coefficient.
                            sigmaSS_TT = @(th) (pi/8)*zetaS^2*delta_rr^2 ...
                                .*(1+cos(th).^2) ...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th)))) * (2*pi*obj.Frequency);

                            sigmaSS_GG = @(th) (pi/8)*zetaS^2*delta_mm^2 ...
                                .*(4*cos(th).^4 - 3*cos(th).^2 + 1) ...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th)))) * (2*pi*obj.Frequency);

                            sigmaSS_GT = @(th) (pi/8)*zetaS^2*rho_mr*delta_mm*delta_rr ...
                                .*(4*cos(th).^3) ...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th)))) * (2*pi*obj.Frequency);

                            sigmaSS = @(th) sigmaSS_TT(th) + sigmaSS_GG(th) + sigmaSS_GT(th);

                            obj.sigma = {sigmaPP, sigmaPS; ...
                                sigmaSP, sigmaSS};
                        case 3
                            %elastic
                            K = obj.vp/obj.vs;
                            zetaP = (2*pi*obj.Frequency)/obj.vp*obj.CorrelationLength;
                            zetaS = K*zetaP;
                            delta_ll = obj.coefficients_of_variation(1); % squared coefficient of variation of lambda
                            delta_mm = obj.coefficients_of_variation(2); % squared coefficient of variation of mu
                            delta_rr = obj.coefficients_of_variation(3); % squared coefficient of variation of density
                            rho_lm = obj.correlation_coefficients(1); % correlation coefficient between lambda and mu
                            rho_lr = obj.correlation_coefficients(2); % correlation coefficient between lambda and density
                            rho_mr = obj.correlation_coefficients(3); % correlation coefficient between mu and density

                            % [Ryzhik et al; Eq. (1.3)] and [Turner, 1998; Eq. (3)]
                            sigmaPP = @(th) (pi/2)*zetaP^3* ...
                                ( (1-2/K^2)^2*delta_ll^2 + 4*(1/K^2-2/K^4)*rho_lm*delta_ll*delta_mm*cos(th).^2 ...
                                + (4/K^4)*delta_mm^2*cos(th).^4 + delta_rr^2*cos(th).^2 ...
                                + 2*(1-2/K^2)*rho_lr*delta_ll*delta_rr*cos(th) ...
                                + (4/K^2)*rho_mr*delta_mm*delta_rr*cos(th).^3 ) ...
                                .*obj.Phi(zetaP.*sqrt(2*(1-cos(th))))*(2*pi*obj.Frequency);

                            %  Ryzhik et al; Eqs. (4.56), (1.20), (1.22)
                            sigmaPS = @(th) (pi/2)*K*zetaP^3* ...
                                ( K^2*delta_rr^2 + 4*delta_mm^2*cos(th).^2 + 4*K*rho_mr*delta_mm*delta_rr*cos(th) )...
                                .*(1-cos(th).^2).*obj.Phi(zetaP.*sqrt(1+K^2-2*K*cos(th)))*(2*pi*obj.Frequency);

                            %  Ryzhik et al; Eqs. (4.56), (1.20), (1.21)
                            sigmaSP = @(th) (pi/4/K^3)*zetaS^3* ...
                                ( delta_rr^2 + (4/K^2)*delta_mm^2*cos(th).^2 + (4/K)*rho_mr*delta_mm*delta_rr*cos(th) ) ...
                                .*(1-cos(th).^2).*obj.Phi(zetaS.*sqrt(1+1/K^2-2/K*cos(th)))*(2*pi*obj.Frequency);

                            % Ryzhik et al; Eq. (4.54)
                            sigmaSS_TT = @(th) (pi/4)*zetaS^3*delta_rr^2*(1+cos(th).^2)...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th))))*(2*pi*obj.Frequency);
                            sigmaSS_GG = @(th) (pi/4)*zetaS^3*delta_mm^2*(4*cos(th).^4-3*cos(th).^2+1)...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th))))*(2*pi*obj.Frequency);
                            sigmaSS_GT = @(th) (pi/4)*zetaS^3*rho_mr*delta_mm*delta_rr*(4*cos(th).^3)...
                                .*obj.Phi(zetaS.*sqrt(2*(1-cos(th))))*(2*pi*obj.Frequency);


                            sigmaSS = @(th) sigmaSS_TT(th) + sigmaSS_GG(th) + sigmaSS_GT(th);

                            obj.sigma = {sigmaPP,sigmaPS; ...
                                sigmaSP,sigmaSS};
                    end
            end
        end
        function calcSigmaAnisotropic(obj)
            %% calcSigmaAnisotropic
            % Compute the normalized differential scattering cross-sections
            % based on statistics of the fluctuating material parameters
            % of the wave equation for anisotropic PSDF
            %
            % Syntax:
            %   newobj = calcSigmaAnisotropic (  );
            %
            % Inputs:
            %
            % Output:

            % The formulas are based on
            % L. Rhyzik, G. Papanicolaou, J. B. Keller. Transport equations for elastic
            % and other waves in random media. Wave Motion 24, pp. 327-370, 1996.
            % doi: 10.1016/S0165-2125(96)00021-2


            error('Not implemented')
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

            % In Ryzhik et al, the PSDFs are those of the fractional parts (normalized)
            % of the corresponding random fields. The corr_matrix thus contains the
            % coefficients of variation in diagonal and correlation coefficients in off-diagonal
            % instead of the standard deviation and covariance, respectively.

            % Following Khazaie et al 2017, PSDFs are supposed as factorizable.
            % Extending to the non-factorizable case is straightforward.

            % constants

            if any(abs(obj.correlation_coefficients)>1)
                error('Absolute values of correlation coefficients should be less than 1!')
            end

            % warning: these are for 3D: formulas should depend on dimensionality
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

            % function to evaluate PSDF:
            % Some notes :
            % The correlation functions are normalized in 1D, 2D, 3D as follows
            % 1D : lc = 2 * integral of R(x)dx from 0 to inf
            % 2D : lc^2 = 2 * integral of x*R(x)dx from 0 to inf
            % 3D : lc^3 = 3 * integral of x^2*R(x)dx from 0 to inf

            % The power spectral density functions are calculated via an n-D
            % Fourier Transform integral, based on the following convention
            % \Phi(k) = (1/2/pi)^d * integral of exp(-i*k*x)*R(x)dx on R^d
            % This integral in 1D, 2D and 3D can be simplified as:
            % 1D : 2/(2*pi) * integral of cos(k*x)*R(x)dx from 0 to inf
            % 2D : 2*pi//(2*pi)^2 * integral of x*J_0(k*x)*R(x)dx from 0 to inf
            % 3D : 4*pi/(2*pi)^3 * integral of x^2*sinc(k*x)*R(x)dx from 0 to inf
        end
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
                obj.R = @(z) exp(-2*sqrt(2)*z);
                obj.Phi = @(z) 1./(16*pi*(1+(z/(2*sqrt(2))).^2).^(1.5));
            elseif obj.d == 3
                obj.R = @(z) exp(-2*6^(1/3)*z);
                obj.Phi = @(z) 1./(48*pi^2*(1+(z/(2*6^(1/3))).^2).^2);
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
            obj.SpectralLaw = 'PowerLaw';
            % out = @(z) 1./(pi^4)*exp(-2*z/pi);
            % obj.Phi = out;
            % obj.R = @(z) 1./(1+(pi^2*z.^2/4))^2;
            if obj.d == 1
                obj.R = @(z) 1./(1+((2/pi)*z)^2)^2;
                obj.Phi = @(z) (pi^3/8)*(z+2/pi).*exp(-pi/2*z);
            elseif obj.d == 2
                obj.R = @(z) 1./(1+(2*z)^2)^2;
                obj.Phi = @(z) (pi/8)*z.*besselk(1,z/2);
            elseif obj.d == 3
                obj.R = @(z) 1./(1+((6*pi)^(1/3)*z)^2)^2;
                obj.Phi = @(z) (pi/6)*exp(-z/((6*pi)^(1/3)));
            else
                error('incorrect dimension!')
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
                obj.Phi = @(z) 1/(2*pi).*exp(-z.^2/4/pi);
            elseif obj.d == 2
                obj.R = @(z) exp(-4*z.^2);
                obj.Phi = @(z) 1./(16*pi)*exp(-z.^2/16);
            elseif obj.d == 3
                obj.R = @(z) exp(-(36*pi)^(1/3)*z.^2);
                obj.Phi = @(z) 1./(48*pi^2)*exp(-z.^2/4/(36*pi)^(2/3));
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
            out = @(z) (3/8/pi^4)*(1-z/2/pi).*obj.heaviside(2*pi-z);
            obj.Phi = out;
            obj.R = @(z)(12*(2-2*cos(2*pi*z)-(2*pi*z).*sin(2*pi*z)))./(2*pi*z).^4;
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
            obj.SpectralLaw = 'LowPass';
            out = @(z) (2/9/pi^4)*obj.heaviside(3*pi/2-z);
            obj.Phi = out;
            obj.R = @(z) (3*(sin(3*pi*z/2)-3*pi*z/2.*cos(3*pi*z/2)))./(3*pi*z/2).^3;
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


            %are the characteristic scales of the medium along the
            % 3-dimensions and $k_x$, $k_y$ and $k_z$
            % are the wavenumber components

            % $K_\nu$ is the modified Bessel function of order $\nu $,
            % where $0.0<\nu<1.0$ is the Hurst number (Mandelbrot, 1985,1983).
            % The fractal dimension of a stochastic field characterized
            % by a von Karman autocorrelation is given by:
            % \begin{displaymath}
            % D=E+1-\nu \
            % \end{displaymath}	(6)
            %
            % where $E$ is the Euclidean dimension i.e., $E=3$
            % for the three-dimensional problem.
            %
            obj.SpectralParam = struct('correlationLength',lc,'nu',nu);
            obj.SpectralLaw = 'VonKarman';
            obj.CorrelationLength = lc;

            if obj.d == 1
                obj.R = @(z) @(z) 2^(1-nu)/gamma(nu) * ...
                    (2*sqrt(pi)*gamma(nu+0.5)/gamma(nu)*z).^nu ...
                    .* besselk(nu,2*sqrt(pi)*gamma(nu+0.5)/gamma(nu)*z);
                obj.Phi = @(z) 1 ./ (1 + (z/(2*sqrt(pi)*gamma(nu+0.5)/gamma(nu))).^2).^(nu+0.5);
            elseif obj.d == 2
                obj.R = @(z) @(z) 2^(1-nu)/gamma(nu) * ...
                    (4*sqrt(nu)*z).^nu .* besselk(nu,4*sqrt(nu)*z);
                obj.Phi = @(z) (pi/4) * 1 ./ (1 + (z/(4*sqrt(nu))).^2).^(nu+0.5);
            elseif obj.d == 3
                obj.R = @(z) @(z) 2^(1-nu)/gamma(nu) * ...
                    ((48*sqrt(pi)*gamma(nu+1.5)/gamma(nu))^(1/3)*z).^nu .* ...
                    besselk(nu,((48*sqrt(pi)*gamma(nu+1.5)/gamma(nu))^(1/3))*z);
                obj.Phi = @(z) (pi/46) * 1 ./ (1 + (z/((48*sqrt(pi)*gamma(nu+1.5)/gamma(nu))^(1/3))).^2).^(nu+1.5);
            else
                error('incorrect dimension!')
            end

            out = obj.Phi;

            % obj.R = @(z) 2^(1-nu)/gamma(nu)* (c*z).^nu .* besselk(nu,c*z);
            % obj.Phi = out;

            % error('Not fully implemented')
            % r = sqrt(x^2/ax^2+y^2/ay^2+z^2/az^2);
            % k=sqrt(kx^2*ax^2+ky^2*ay^2+kz^2*az^2);
            % bessel = besseli(nu,r);
            % C = 4*pi*nu*H^2*r^nu*bessel/bessel(1);
            % S = (4*pi*nu*H^2/bessel(1))*(ax^2+ay^2+az^2)/((1+k^2^(nu+1.5)));
            % out = @(z) 1;
            % obj.Phi = out;
            % obj.R = out;

        end
        function out = MonoDisperseSphere(obj, eta, D)
            %% MonoDisperseSphere
            % Compute the normalized power spectral density function for a
            % monodisperse hard-sphere suspension (3D) and, additionally,
            % solve the 2D Ornstein-Zernike equation with the Percus-Yevick
            % closure for an equivalent hard-disk system to obtain the full
            % set of microstructural descriptors:
            %   g(r)   - pair correlation function (radial distribution function)
            %   S(k)   - static structure factor
            %   S2(r)  - two-point probability function
            %   chi(r) - autocovariance  chi(r) = S2(r) - phi^2
            %
            % The 3D part follows:
            %   Torquato S. and Stell G. (1985). Microstructure of two-phase
            %   random media. V. The n-point matrix probability functions for
            %   impenetrable spheres, J. Chem. Phys., 82(2), pp. 980-987.
            %
            % The 2D OZ+PY part mirrors PY2D_disk.m:
            %   Adda-Bedia, Katzav & Vella, J. Chem. Phys. 128, 184508 (2008)
            %   Torquato, Random Heterogeneous Materials, Springer (2002)
            %
            % Syntax:
            %   out = MonoDisperseSphere(obj, eta, D)
            %
            % Inputs:
            %   eta : packing fraction (area fraction in 2D / volume fraction in 3D)
            %   D   : sphere/disk diameter (same units as desired length scale)
            %
            % Output:
            %   out : function handle for the PSDF  Phi(k)  (3D, normalized)
            %
            % Side effects:
            %   obj.R, obj.Phi, obj.k, obj.CorrelationLength are set (3D part).
            %   A 2x2 figure panel is produced with g(r), S(k), S2(r), chi(r)
            %   computed from the 2D PY solver (analogous to PY2D_disk/plot_results).

            % =================================================================
            % PART 1 — 3D monodisperse sphere PSDF  (original implementation)
            % =================================================================

            switch obj.d
                case 1
                    error('Power spectral density for monodisperse sphere not supported in 1D')
                case 2
                    % =================================================================
                    % PART 2 — 2D Percus-Yevick OZ solver  (from PY2D_disk.m)
                    %          Produces g(r), S(k), S2(r), chi(r) for the equivalent
                    %          hard-disk system and plots all four panels.
                    % =================================================================
                    fprintf('\n--- 2D PY OZ solver (MonoDisperseSphere) ---\n');
                    fprintf('    phi = %.4f,  D = %.6g\n', eta, D);

                    % --- Grid parameters ---
                    Nr      = 4096;
                    r_max   = 20.0 * D;
                    dr      = r_max / Nr;
                    r_py    = ((1:Nr) - 0.5)' * dr;   % cell-centred, avoids r = 0

                    dk_py   = pi / r_max;
                    k_max   = pi / dr;
                    Nk_py   = floor(k_max / dk_py);
                    k_py    = ((1:Nk_py) - 0.5)' * dk_py;

                    % --- Picard iteration parameters ---
                    rho_2D   = eta / (pi*(D/2)^2);   % 2D number density
                    inside   = r_py < D;
                    outside  = ~inside;
                    max_iter = 5000;
                    tol      = 1e-12;
                    alpha    = 0.4;

                    % Hankel transform helpers (defined as nested functions below)
                    hankel_fwd = @(f) MaterialClass.hankel0_fwd(f, r_py, k_py, dr);
                    hankel_inv = @(f) MaterialClass.hankel0_inv(f, k_py, r_py, dk_py);

                    % Initial guess
                    gamma_r = zeros(Nr, 1);
                    c_r     = zeros(Nr, 1);
                    c_r(inside) = -1.0;

                    converged = false;
                    for iter = 1:max_iter
                        c_k         = hankel_fwd(c_r);
                        h_k         = c_k ./ (1.0 - rho_2D * c_k);
                        h_r         = hankel_inv(h_k);
                        gamma_r_new = h_r - c_r;
                        err         = max(abs(gamma_r_new - gamma_r));
                        gamma_r     = (1 - alpha)*gamma_r + alpha*gamma_r_new;
                        c_r(inside)  = -1.0 - gamma_r(inside);
                        c_r(outside) =  0.0;
                        if err < tol
                            fprintf('    PY converged at iter %d, err = %.3e\n', iter, err);
                            converged = true;
                            break;
                        end
                    end
                    if ~converged
                        fprintf('    WARNING: PY did not fully converge (err=%.3e).\n', err);
                    end

                    % Final correlation functions
                    c_k  = hankel_fwd(c_r);
                    h_k  = c_k ./ (1.0 - rho_2D * c_k);
                    h_r  = hankel_inv(h_k);
                    g_r  = 1.0 + h_r;
                    g_r(inside) = 0.0;              % enforce hard core
                    S_k  = 1.0 + rho_2D * h_k;     % static structure factor

                    % Two-point probability function S2(r) via spectral density route
                    % (Torquato & Lu approach):  chi_hat(k) = rho * |m_hat(k)|^2 * S(k)
                    R_disk  = D / 2;
                    kR      = k_py * R_disk;
                    m_hat   = zeros(Nk_py, 1);
                    idx0    = kR < 1e-8;
                    m_hat(idx0)  = pi * R_disk^2;
                    m_hat(~idx0) = 2*pi*R_disk^2 * besselj(1, kR(~idx0)) ./ kR(~idx0);

                    chi_hat = rho_2D * abs(m_hat).^2 .* S_k;
                    chi_r   = hankel_inv(chi_hat);
                    S2_r    = max(0, min(1, eta^2 + chi_r));

                    % Equation of state printout
                    [~, idx_c] = min(abs(r_py - D));
                    g_contact  = g_r(idx_c);
                    Z_virial   = 1 + (rho_2D/2) * pi * D^2 * g_contact;
                    fprintf('    g(D+) = %.6f,  Z_virial = %.6f\n', g_contact, Z_virial);
                    fprintf('    S(k->0) = %.6f  =>  Z_comp ~ %.6f\n', S_k(1), 1/S_k(1));

                    % --- 2x2 figure: same layout as plot_results in PY2D_disk.m ---
                    r_norm_py  = r_py  / D;
                    k_norm_py  = k_py  * D;
                    r_plot_max = 6.0;
                    k_plot_max = 20.0;
                    idx_r      = r_norm_py <= r_plot_max;
                    idx_k      = k_norm_py <= k_plot_max;

                    fig = figure('Name','PY 2D Hard Disks — MonoDisperseSphere', ...
                        'Color','w', 'Position',[100 80 1200 900]);
                    tl = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
                    title(tl, sprintf('2D Percus-Yevick — Hard Disks  (\\phi = %.2f,  D = %.4g)', ...
                        eta, D), 'FontSize', 14, 'FontWeight','bold');

                    % g(r)
                    nexttile;
                    plot(r_norm_py(idx_r), g_r(idx_r), 'b-', 'LineWidth', 1.8);
                    hold on;
                    yline(1, 'k--', 'LineWidth', 0.8);
                    xline(1, 'r:', 'D', 'LineWidth', 0.8, 'LabelVerticalAlignment','bottom');
                    xlabel('r / D');  ylabel('g(r)');
                    title('Pair Correlation Function g(r)');
                    xlim([0 r_plot_max]);
                    ylim([0 max(g_r(idx_r))*1.05 + 0.1]);
                    grid on;  box on;

                    % S(k)
                    nexttile;
                    plot(k_norm_py(idx_k), S_k(idx_k), 'r-', 'LineWidth', 1.8);
                    hold on;
                    yline(1, 'k--', 'LineWidth', 0.8);
                    xlabel('k D');  ylabel('S(k)');
                    title('Static Structure Factor S(k)');
                    xlim([0 k_plot_max]);
                    ylim([0 max(S_k(idx_k))*1.05 + 0.1]);
                    grid on;  box on;

                    % S2(r)
                    nexttile;
                    plot(r_norm_py(idx_r), S2_r(idx_r), 'm-', 'LineWidth', 1.8);
                    hold on;
                    yline(eta^2, 'k--', 'LineWidth', 0.8, 'Label','\phi^2', ...
                        'LabelVerticalAlignment','bottom');
                    yline(eta,   'k:',  'LineWidth', 0.8, 'Label','\phi');
                    xline(1, 'r:', 'D', 'LineWidth', 0.8, 'LabelVerticalAlignment','bottom');
                    xlabel('r / D');  ylabel('S_2(r)');
                    title('Two-Point Correlation S_2(r)');
                    xlim([0 r_plot_max]);
                    ylim([0 eta + 0.05]);
                    grid on;  box on;

                    % chi(r)
                    nexttile;
                    plot(r_norm_py(idx_r), chi_r(idx_r), 'g-', 'LineWidth', 1.8);
                    hold on;
                    yline(0, 'k--', 'LineWidth', 0.8);
                    xline(1, 'r:', 'D', 'LineWidth', 0.8, 'LabelVerticalAlignment','bottom');
                    xlabel('r / D');  ylabel('\chi(r) = S_2(r) - \phi^2');
                    title('Autocovariance \chi(r)');
                    xlim([0 r_plot_max]);
                    grid on;  box on;
                case 3
                    % discretization in Fourier space (k normalised to sphere radius = D/2 = 1)
                    obj.k = linspace(0, 6/D, 4094);
                    % discretization in real space
                    raux = linspace(0, 5*D, 4096);

                    % number density of spheres  rhoS = phi / V_sphere
                    rhoS = 3*eta / (4*pi);   % with radius normalised to 1 (r_norm = 2r/D)

                    % normalise radii to sphere radius = 1
                    raux = 2*raux / D;

                    % Union volume of two unit-radius spheres at separation r (normalised)
                    V2 = 8*pi/3 * ones(size(raux));
                    V2(raux < 2) = 4*pi/3 * (1 + 3/4*raux(raux<2) - raux(raux<2).^3/16);

                    % PY direct correlation function c(k) — analytic solution for 3D hard spheres
                    l1 = (1 + 2*eta)^2 / (1 - eta)^4;
                    l2 = -(1 + eta/2)^2 / (1 - eta)^4;

                    c = -4*pi ./ (obj.k.^3) .* ( ...
                        l1*(sin(2*obj.k) - 2*obj.k.*cos(2*obj.k)) + ...
                        3*eta*l2./obj.k .* (4*obj.k.*sin(2*obj.k) + (2 - 4*obj.k.^2).*cos(2*obj.k) - 2) + ...
                        eta*l1./(2*obj.k.^3) .* ((6*obj.k.^2 - 3 - 2*obj.k.^4).*cos(2*obj.k) + ...
                        (4*obj.k.^3 - 6*obj.k).*sin(2*obj.k) + 3) );
                    c(1) = -8*pi/3 * ((4 + eta)*l1 + 18*eta*l2);

                    % OZ relation: h_hat = c_hat / (1 - rhoS*c_hat)
                    h3D = c ./ (1 - rhoS*c);

                    % Fourier transform of the sphere indicator function m_hat(k)
                    m3D = 4*pi * ((sin(obj.k)./obj.k - cos(obj.k)) ./ obj.k.^2);
                    m3D(1) = 4*pi/3;

                    % 2-point matrix probability function M(r) via inverse sine transform
                    M = zeros(size(raux));
                    M(1) = -(eta/rhoS)^2;
                    for i1 = 2:length(raux)
                        M(i1) = 1/2/pi^2/raux(i1) * trapz(obj.k, h3D.*m3D.^2.*obj.k.*sin(obj.k*raux(i1)));
                    end

                    S2_3D = 1 - rhoS*V2 + rhoS^2*M + eta^2;

                    % Particle autocorrelation (normalised)
                    Racf = (S2_3D - (1-eta)^2) / eta / (1-eta);

                    % Convert back to physical units and apply Hann window at tail
                    raux = D*raux/2;
                    right_window_length = 200;
                    hann_win = hann(2*right_window_length);
                    Racf(end-right_window_length+1:end) = ...
                        Racf(end-right_window_length+1:end) .* hann_win(right_window_length+1:end)';

                    % Correlation length (3D definition: lc^3 = 3 * int r^2 R(r) dr)
                    obj.R = @(z) interp1(raux, Racf, z, 'makima', 0);
                    LL = 2*integral(obj.R, 0, inf);
                    obj.CorrelationLength = LL;
                    raux_norm = raux / LL;
                    obj.R = @(z) interp1(raux_norm, Racf, z, 'makima', 0);

                    % PSDF via 3D Fourier-Bessel transform
                    obj.k = obj.k * obj.CorrelationLength;
                    phi_k = zeros(1, length(obj.k));
                    for i1 = 1:length(obj.k)
                        phi_k(i1) = 1/2/pi^2 * trapz(raux_norm, ...
                            sinc(raux_norm * obj.k(i1)) .* raux_norm.^2 .* Racf);
                    end
                    P = abs(phi_k .* conj(phi_k));
                    obj.Phi = @(z) interp1(obj.k, P, z, 'makima', 0);
                    out = obj.Phi;
            end
        end
        %% function to evaluate PSDF: an image
        function out = GetPSDFromImage(obj, Im, dx, dy, dz)
            %% GetPSDFromImage
            % Compute the normalised power spectral density function from a
            % 2D grayscale image or a 3D binary volume (logical/uint8/double).
            %
            % The method replaces the previous 1-D xcorr approach with the
            % exact 3D FFT-based pipeline translated from Center2S2func.py:
            %   1. Convert input to a binary (0/1) indicator field.
            %   2. Compute S2(r) via the Wiener-Khinchin theorem:
            %         S2_map = IFFT( |FFT(I)|^2 ) / N_voxels
            %   3. Radially average S2_map to obtain the isotropic S2(r).
            %   4. Derive the normalised autocorrelation (ACF):
            %         R(r) = [S2(r) - phi^2] / [phi*(1-phi)]
            %   5. Compute the correlation length Lc from R(r).
            %   6. Compute the 3D PSDF via Fourier-Bessel quadrature and
            %      store obj.Phi, obj.R, obj.k, obj.CorrelationLength.
            %
            % The method also calls CalcS2Correlation (static) so that the
            % full S2 pipeline is always available as a standalone utility.
            %
            % Syntax:
            %   out = GetPSDFromImage(obj, Im)
            %   out = GetPSDFromImage(obj, Im, dx)
            %   out = GetPSDFromImage(obj, Im, dx, dy)
            %   out = GetPSDFromImage(obj, Im, dx, dy, dz)
            %
            % Inputs:
            %   Im  : 2D or 3D numeric array (grayscale image or binary volume).
            %         Values > 0.5 are treated as the solid phase (indicator = 1).
            %   dx  : voxel/pixel size along x  [same unit as desired Lc] (default 1)
            %   dy  : voxel/pixel size along y  (default dx)
            %   dz  : voxel/pixel size along z  (default dx, ignored for 2D)
            %
            % Output:
            %   out : function handle  Phi(k)  — normalised 3D PSDF
            %
            % Side effects (properties set on obj):
            %   obj.Phi              — PSDF function handle
            %   obj.R                — normalised ACF function handle
            %   obj.k                — wavenumber vector (normalised by Lc)
            %   obj.r                — radial distance vector (normalised by Lc)
            %   obj.CorrelationLength— correlation length Lc  [same units as dx]
            %
            % See also: MaterialClass.CalcS2Correlation, MaterialClass.VoxelizeDomain

            % --- default pixel sizes ---
            if nargin < 3 || isempty(dx),  dx = 1;  end
            if nargin < 4 || isempty(dy),  dy = dx; end
            if nargin < 5 || isempty(dz),  dz = dx; end

            % --- convert to binary double indicator field ---
            if islogical(Im)
                I = double(Im);
            elseif ndims(Im) == 3 && size(Im,3) == 3
                % RGB image — convert to grayscale first
                I = double(rgb2gray(Im)) / 255;
                I = double(I > 0.5);
            else
                I = double(Im);
                if max(I(:)) > 1,  I = I / max(I(:));  end
                I = double(I > 0.5);
            end

            % domain size in physical units
            [nx, ny, nz] = size(I);
            if nz == 1
                % 2D image: treat as single-voxel-thick volume
                L = [nx*dx, ny*dy, dz];
                resolucao = min([dx dy]);
            else
                L = [nx*dx, ny*dy, nz*dz];
                resolucao = min([dx dy dz]);
            end

            % --- S2 via FFT (CalcS2Correlation) ---
            [r_axis, S2_radial, ~, phi_vol, ~] = ...
                MaterialClass.CalcS2Correlation(I, resolucao, L);

            % --- normalised autocorrelation R(r) = [S2(r)-phi^2]/[phi*(1-phi)] ---
            % Avoid division by zero for degenerate (all-0 or all-1) images
            denom = phi_vol * (1 - phi_vol);
            if denom < 1e-12
                warning('MaterialClass:GetPSDFromImage', ...
                    'Volume fraction is %.4f — image may be degenerate.', phi_vol);
                denom = 1;
            end
            Racf = (S2_radial - phi_vol^2) / denom;
            Racf = Racf(:);          % ensure column vector
            r_phys = r_axis(:);     % physical radial axis [same unit as dx]

            % Apply Hann window to suppress tail oscillations
            right_win = min(200, floor(numel(Racf)/4));
            win = hann(2*right_win);
            Racf(end-right_win+1:end) = ...
                Racf(end-right_win+1:end) .* win(right_win+1:end);

            % --- correlation length (3D definition: lc^3 = 3*int r^2 R dr) ---
            % Use the simpler integral form: Lc = 2 * int_0^inf R(r) dr
            R_interp = @(z) interp1(r_phys, Racf, z, 'makima', 0);
            Lc = 2 * integral(R_interp, 0, r_phys(end));
            if Lc <= 0
                warning('MaterialClass:GetPSDFromImage', ...
                    'Correlation length came out non-positive (%.4g). Check image.', Lc);
                Lc = r_phys(end) / 10;
            end
            obj.CorrelationLength = Lc;

            % normalise r by Lc
            r_norm = r_phys / Lc;
            obj.r  = r_norm;
            obj.R  = @(z) interp1(r_norm, Racf, z, 'makima', 0);

            % --- 3D PSDF via Fourier-Bessel (sinc) quadrature ---
            %   Phi(k) = 4*pi/(2*pi)^3 * int_0^inf r^2 sinc(k*r) R(r) dr
            %   with r and k both normalised by Lc
            Nk_psd = 512;
            k_max_psd = pi / mean(diff(r_norm));  % Nyquist in normalised units
            k_max_psd = min(k_max_psd, 6.0);       % cap at 6/Lc (physically meaningful)
            obj.k = linspace(0, k_max_psd, Nk_psd);

            psd_vals = zeros(1, Nk_psd);
            for ik = 1:Nk_psd
                if obj.k(ik) < 1e-8
                    integrand = r_norm.^2 .* Racf;
                else
                    kr = obj.k(ik) * r_norm;
                    j0 = sin(kr) ./ kr;
                    j0(kr < 1e-12) = 1;
                    integrand = r_norm.^2 .* Racf .* j0;
                end
                psd_vals(ik) = (4*pi / (2*pi)^3) * trapz(r_norm, integrand);
            end
            psd_vals = max(psd_vals, 0);   % enforce non-negativity

            obj.Phi = @(z) interp1(obj.k, psd_vals, z, 'makima', 0);
            out = obj.Phi;

            % --- summary figure: S2(r), R(r), PSDF ---
            r_plot  = r_norm(r_norm <= 6);
            S2_plot = S2_radial(r_norm <= 6);
            R_plot  = Racf(r_norm <= 6);

            fig = figure('Name','GetPSDFromImage — Microstructure Descriptors', ...
                'Color','w', 'Position',[120 80 1200 400]);
            tl = tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');
            title(tl, sprintf('Image PSDF  (\\phi = %.4f,  L_c = %.4g)', phi_vol, Lc), ...
                'FontSize', 13, 'FontWeight', 'bold');

            nexttile;
            plot(r_plot, S2_plot, 'm-', 'LineWidth', 1.8);
            hold on;
            yline(phi_vol^2, 'k--', 'LineWidth', 0.8, 'Label', '\phi^2', ...
                'LabelVerticalAlignment','bottom');
            yline(phi_vol,   'k:',  'LineWidth', 0.8, 'Label', '\phi');
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
            plot(obj.k, psd_vals, 'r-', 'LineWidth', 1.8);
            xlabel('k \cdot L_c');  ylabel('\Phi(k)');
            title('Power Spectral Density \Phi(k)');
            xlim([0 min(6, obj.k(end))]);  grid on;  box on;
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
            
            % Use a logarithmic grid for integration to capture decay
            r = logspace(-4, 3, 4096);  % covers range from very small to large
            R_vals = obj.R(r);
            
            if obj.d == 1
                % Lc = 2 * int_0^inf R(r) dr
                Lc = 2 * trapz(r, R_vals);
            elseif obj.d == 2
                % Lc^2 = 2 * int_0^inf r*R(r) dr
                I = trapz(r, r .* R_vals);
                Lc = sqrt(2 * I);
            elseif obj.d == 3
                % Lc^3 = 3 * int_0^inf r^2*R(r) dr
                I = trapz(r, r.^2 .* R_vals);
                Lc = (3 * I)^(1/3);
            else
                error('MaterialClass:InvalidDimension', 'Dimension must be 1, 2, or 3');
            end
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
            r = logspace(-4, 3, 4096);  % real space grid
            R_vals = obj.R(r);
            
            % Wavenumber grid
            k = logspace(-3, 2, 512);  % wavenumber grid
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
            obj.r = r;
            obj.k = k;
            obj.Phi = @(z) interp1(k, PSD, z, 'pchip', 0);
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
                    % Pick arbitrary tangent [1 0 0] projected on plane
                    nUnit = n(:).';
                    t = repmat([1 0 0], sum(nearNormal), 1);
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
            % Ideally, call this once and store 'out' in the object properties
            [out, angles] = MaterialClass.Zoeppritz(obj);

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
        %% OTHER
        function h = heaviside(h)
            h(h>0) = 1;
            h(h==0) = 1/2;
            h(h<0) = 0;
        end
    end
    methods(Static)
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
            %  mat: MaterialClass object or a vector of MaterialClass
            %  object.
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
                % here we should make a combination of all the material
                % interfaces and evaluate...
                error("Not implemented")
                % angles
                j1_deg = linspace(0,90,181);

                if(mat.acoustics)
                    error("This dont work for acoustic material")
                end

                vp1  = mat(1).vp;
                vs1  = mat(1).vs;
                rho1 = mat(1).rho;

                vp2  = mat(2).vp;
                vs2  = mat(2).vs;
                rho2 = mat(2).rho;

                if vs2 == 0
                    out = MaterialClass.ZoeppritzFluid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
                else
                    out = MaterialClass.ZoeppritzSolid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
                end

                o.Rsh   = @(z)interp1(j1_deg,out.E_Rsh,z);
                o.Tsh   = @(z)interp1(j1_deg,out.E_Tsh,z);

                o.Rpp   = @(z)interp1(j1_deg,out.E_Rpp ,z);
                o.Rpsv  = @(z)interp1(j1_deg,out.E_Rpsv,z);
                o.Tpp   = @(z)interp1(j1_deg,out.E_Tpp ,z);
                o.Tpsv  = @(z)interp1(j1_deg,out.E_Tpsv,z);

                o.Rsvp   = @(z)interp1(j1_deg,out.E_Rsp  ,z);
                o.Rsvsv = @(z)interp1(j1_deg,out.E_Rsvsv,z);
                o.Tsvp   = @(z)interp1(j1_deg,out.E_Tsp  ,z);
                o.Tsvsv = @(z)interp1(j1_deg,out.E_Tsvsv,z);
            end
        end
        out = ZoeppritzFluid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
        out = ZoeppritzSolid(j1_deg,vp1,vs1,rho1,vp2,vs2,rho2);
        function Microestrutura = VoxelizeDomain(coordinates, D, L, resolucao)
            %% VoxelizeDomain
            % Discretise a periodic 3D domain into a binary voxel grid and
            % paint spheres centred at the given coordinates.
            %
            % This is a direct MATLAB translation of voxelize_domain() from
            % Center2S2func.py, keeping the same periodic-boundary logic.
            %
            % Syntax:
            %   Microestrutura = MaterialClass.VoxelizeDomain(coordinates, D, L, resolucao)
            %
            % Inputs:
            %   coordinates : (N x 3) array of sphere centre positions [x y z]
            %                 in the same units as D and L.
            %   D           : sphere diameter.
            %   L           : 1x3 domain size  [Lx  Ly  Lz].
            %   resolucao   : voxel edge length (isotropic).
            %
            % Output:
            %   Microestrutura : logical array of size [nx ny nz].
            %                    true  = solid (sphere) phase.
            %                    false = matrix phase.
            %
            % Example:
            %   % 200 random spheres, D=4, box 100x100x100, voxel=1
            %   ctrs = rand(200,3) .* [100 100 100];
            %   M = MaterialClass.VoxelizeDomain(ctrs, 4, [100 100 100], 1);
            %   phi_vox = mean(M(:))   % should be close to 200*pi/6*4^3/1e6

            fprintf('Voxelizando o dominio...\n');

            % domain discretisation
            nx = ceil(L(1) / resolucao);
            ny = ceil(L(2) / resolucao);
            nz = ceil(L(3) / resolucao);

            % initialise as matrix phase (false)
            Microestrutura = false(nx, ny, nz);

            R = D / 2;   % sphere radius

            % voxel-centre coordinate vectors
            xv = linspace(0, L(1), nx);
            yv = linspace(0, L(2), ny);
            zv = linspace(0, L(3), nz);

            n_spheres = size(coordinates, 1);
            report_step = max(1, floor(n_spheres / 10));

            for i = 1:n_spheres
                % progress report every ~10 %
                if n_spheres > 100 && mod(i, report_step) == 0
                    fprintf('  Progresso: %.1f%%\n', 100*i/n_spheres);
                end

                cx = coordinates(i,1);
                cy = coordinates(i,2);
                cz = coordinates(i,3);

                % indices within bounding box (flat distance <= R)
                ix = find(abs(xv - cx) <= R);
                iy = find(abs(yv - cy) <= R);
                iz = find(abs(zv - cz) <= R);

                % wrap for periodic boundaries (1-based indexing)
                ix = mod(ix - 1, nx) + 1;
                iy = mod(iy - 1, ny) + 1;
                iz = mod(iz - 1, nz) + 1;

                if isempty(ix) || isempty(iy) || isempty(iz)
                    continue;
                end

                % 3D meshgrid of candidate voxel indices
                [IX, IY, IZ] = ndgrid(ix, iy, iz);

                % physical coordinates of candidate voxels
                Xc = xv(IX);
                Yc = yv(IY);
                Zc = zv(IZ);

                % signed distances with periodic correction
                ddx = abs(Xc - cx);   ddx = ddx - (ddx > L(1)/2) * L(1);
                ddy = abs(Yc - cy);   ddy = ddy - (ddy > L(2)/2) * L(2);
                ddz = abs(Zc - cz);   ddz = ddz - (ddz > L(3)/2) * L(3);

                % sphere mask
                mask = ddx.^2 + ddy.^2 + ddz.^2 <= R^2;

                % linear indices for assignment
                lin_idx = sub2ind([nx ny nz], IX(mask), IY(mask), IZ(mask));
                Microestrutura(lin_idx) = true;
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
            % The 3D map is then radially averaged to give the 1-D function S2(r).
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

                % Diffusion coefficient m²/s (Eq. (5.12), Ryzhik et al, 1996)
                mat.Diffusivity = mat.v^2/(double(d)*(mat.Sigma-mat.Sigmapr));

                mat.meanFreeTime = 1/mat.Sigma;
                mat.meanFreePath = mat.v * mat.meanFreeTime;

                mat.transportMeanFreePath = double(d) * mat.Diffusivity / mat.v;
                mat.transportMeanFreeTime = mat.transportMeanFreeTime;

                % anisotropy coefficient (characterizes scattering directionality)
                mat.g = 1 - mat.meanFreePath/mat.transportMeanFreePath;

                % the two lines below are just for homogenization of the propagation
                % code between acoustics and elastics
                mat.vp = mat.v;
                mat.vs = 0;
            else
                mat.Sigma = zeros(2);
                mat.invcdf = cell(2);

                [mat.Sigma(1,1),mat.Sigmapr(1,1),mat.invcdf{1,1}] = MaterialClass.prepareSigmaOne(mat.sigma{1,1},d);
                [mat.Sigma(1,2),mat.Sigmapr(1,2),mat.invcdf{1,2}] = MaterialClass.prepareSigmaOne(mat.sigma{1,2},d);
                [mat.Sigma(2,1),mat.Sigmapr(2,1),mat.invcdf{2,1}] = MaterialClass.prepareSigmaOne(mat.sigma{2,1},d);
                [mat.Sigma(2,2),mat.Sigmapr(2,2),mat.invcdf{2,2}] = MaterialClass.prepareSigmaOne(mat.sigma{2,2},d);
                mat.meanFreeTime = 1./sum(mat.Sigma,2);
                mat.meanFreePath = [mat.vp mat.vs].' .* mat.meanFreeTime;

                K = mat.vp/mat.vs;
                % Transport mean free paths of P & S waves
                tmfp_P = (mat.vp*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2)) + mat.vs*mat.Sigmapr(1,2) )...
                                  /( (mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1))*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2))- mat.Sigmapr(1,2)*mat.Sigmapr(2,1) );
                tmfp_S = (mat.vs*(mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1)) + mat.vp*mat.Sigmapr(2,1) )...
                                  /( (mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1))*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2))- mat.Sigmapr(1,2)*mat.Sigmapr(2,1) );
                mat.transportMeanFreePath = [tmfp_P tmfp_S]';
                mat.transportMeanFreeTime = mat.transportMeanFreePath / [mat.vp mat.vs].';

                % partial diffusion coefficients of P & S waves
                Dp = mat.vp*tmfp_P/d;
                Ds = mat.vs*tmfp_S/d;
                % Diffusion coefficient m²/s (Eqs. (5.42) & (5.46), Ryzhik et al, 1996)
                mat.Diffusivity = double((Dp+2*K^3*Ds)/(1+2*K^3));

                mat.P2P = mat.Sigma(1,1)/sum(mat.Sigma(1,:),2);
                mat.S2S = mat.Sigma(2,2)/sum(mat.Sigma(2,:),2);
            end
        end
        function [Sigma,Sigma_prime,invcdf] = prepareSigmaOne(sigma,d)
            Nth = 1e6;
            xth = linspace(0,pi,Nth);
            if d==2
                Sigma = 2*integral(sigma,0,pi);
                Sigma_prime = 2*integral(@(th)sigma(th).*cos(th),0,pi);
                sigmaNorm = @(th) (2/Sigma)*sigma(th);
            elseif d==3
                Sigma = 2*pi*integral(@(th)sigma(th).*sin(th),0,pi);
                Sigma_prime = 2*pi*integral(@(th)sigma(th).*sin(th).*cos(th),0,pi);
                sigmaNorm = @(th) (2*pi/Sigma)*sin(th).*sigma(th);
            end
            pdf = sigmaNorm(xth);
            cdf = cumsum(pdf)*(pi/Nth);
            ind = find(diff(cdf)>1e-8);
            ind = unique([ind ind+1]);
            invcdf = griddedInterpolant(cdf(ind),xth(ind));
        end
    end
end