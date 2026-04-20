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
        SpectralLaw_def = {'','exp','power_law','gaussian','triangular','low_pass','VonKarman','monodispersesphere','image','Imported'};
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
                if geometry.dimension==1
                    warning(['Total scattering cross-sections for 1D case are to be ' ...
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
            obj.SpectralLaw = 'power_law';
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
                obj.R = @(z) 1./(1+((6*pi)^(1/3)*z).^2).^2;
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
            obj.SpectralLaw = 'low_pass';
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
                obj.R = @(z) 2^(1-nu)/gamma(nu) * ...
                    (2*sqrt(pi)*gamma(nu+0.5)/gamma(nu)*z).^nu ...
                    .* besselk(nu,2*sqrt(pi)*gamma(nu+0.5)/gamma(nu)*z);
                obj.Phi = @(z) 1 ./ (1 + (z./(2*sqrt(pi)*gamma(nu+0.5)/gamma(nu))).^2).^(nu+0.5);
            elseif obj.d == 2
                obj.R = @(z) 2^(1-nu)/gamma(nu) * ...
                    (4*sqrt(nu)*z).^nu .* besselk(nu,4*sqrt(nu)*z);
                obj.Phi = @(z) (pi/4) * 1 ./ (1 + (z./(4*sqrt(nu))).^2).^(nu+0.5);
            elseif obj.d == 3
                obj.R = @(z) 2^(1-nu)/gamma(nu) * ...
                    ((48*sqrt(pi)*gamma(nu+1.5)/gamma(nu))^(1/3)*z).^nu .* ...
                    besselk(nu,((48*sqrt(pi)*gamma(nu+1.5)/gamma(nu))^(1/3))*z);
                obj.Phi = @(z) (pi/46) * 1 ./ (1 + (z./((48*sqrt(pi)*gamma(nu+1.5)/gamma(nu)).^(1/3))).^2).^(nu+1.5);
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

                    
                    Racf = (S2_r-eta^2)./(eta*(1-eta))

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
                    r = 0:D/40:25*D;
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

                    obj.r = r;
                    obj.R = @(z) interp1(r, Racf, z, 'makima', 0);
                    obj.CalcLc;
                    r_norm = r / obj.CorrelationLength;
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
            % 1-D branch: vector input
            % =========================================================
            if isvector(I)
                I       = I(:);
                N       = numel(I);
                phi_vol = mean(I);

                % S2 via FFT (Wiener-Khinchin), periodic domain
                F      = fft(I);
                S2_raw = real(ifft(F .* conj(F))) / N;
                S2_raw = fftshift(S2_raw);          % r=0 at centre

                % take non-negative-r half
                half      = floor(N/2);
                S2_radial = S2_raw(half+1 : end);
                r_axis    = (0 : numel(S2_radial)-1)' * dx;

                % clip to L/2
                keep      = r_axis <= N*dx/2;
                r_axis    = r_axis(keep);
                S2_radial = S2_radial(keep);

                % normalised ACF
                denom = phi_vol * (1 - phi_vol);
                if denom < 1e-12
                    warning('MaterialClass:GetPSDFromImage', ...
                        'Volume fraction is %.4f — image may be degenerate.', phi_vol);
                    denom = 1;
                end
                Racf   = (S2_radial - phi_vol^2) / denom;
                r_phys = r_axis;

                % Hann window
                right_win = min(200, floor(numel(Racf)/4));
                win = hann(2*right_win);
                Racf(end-right_win+1:end) = ...
                    Racf(end-right_win+1:end) .* win(right_win+1:end);

                % correlation length: Lc = 2 * int_0^inf R(r) dr  (1D definition)
                Lc = 2 * trapz(r_phys, max(Racf, 0));
                if Lc <= 0
                    warning('MaterialClass:GetPSDFromImage', ...
                        'Correlation length came out non-positive (%.4g). Check image.', Lc);
                    Lc = r_phys(end) / 10;
                end
                obj.CorrelationLength = Lc;

                r_norm = r_phys / Lc;
                obj.r  = r_norm;
                obj.R  = @(z) interp1(r_norm, Racf, z, 'makima', 0);

                % 1D PSDF: Phi(k) = (1/pi) * int_0^inf R(r) cos(k*r) dr
                Nk_psd    = 512;
                k_max_psd = min(pi / mean(diff(r_norm)), 6.0);
                obj.k     = linspace(0, k_max_psd, Nk_psd);
                psd_vals  = zeros(1, Nk_psd);
                for ik = 1:Nk_psd
                    psd_vals(ik) = (1/pi) * trapz(r_norm, Racf .* cos(obj.k(ik) .* r_norm));
                end
                psd_vals = max(psd_vals, 0);
                obj.Phi  = @(z) interp1(obj.k, psd_vals, z, 'makima', 0);
                out      = obj.Phi;

                r_plot  = r_norm(r_norm <= 6);
                S2_plot = S2_radial(r_norm <= 6);
                R_plot  = Racf(r_norm <= 6);
                MaterialClass.psd_summary_figure('1D', phi_vol, Lc, ...
                    r_plot, S2_plot, R_plot, obj.k, psd_vals, phi_vol);
                return;
            end

            % =========================================================
            % 2-D / 3-D branch
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
            Racf   = (S2_radial - phi_vol^2) / denom;
            Racf   = Racf(:);
            r_phys = r_axis(:);

            % correlation length — set obj.r/obj.R in physical units first so
            % CalcLc can integrate using the correct dimension formula (obj.d)
            obj.r = r_phys;
            obj.R = @(z) interp1(r_phys, Racf, z, 'makima', 0);
            Lc = obj.CalcLc;
            if Lc <= 0
                warning('MaterialClass:GetPSDFromImage', ...
                    'Correlation length came out non-positive (%.4g). Check image.', Lc);
                Lc = r_phys(end) / 10;
            end
            obj.CorrelationLength = Lc;

            % re-normalise r by Lc
            r_norm = r_phys / Lc;
            obj.r  = r_norm;
            obj.R  = @(z) interp1(r_norm, Racf, z, 'makima', 0);

            % PSDF via Fourier-Bessel quadrature
            %   2D: Phi(k) = 1/(2*pi) * int_0^inf r*J0(k*r)*R(r) dr
            %   3D: Phi(k) = 4*pi/(2*pi)^3 * int_0^inf r^2*sinc(k*r)*R(r) dr
            Nk_psd    = 512;
            k_max_psd = min(pi / mean(diff(r_norm)), 6.0);
            obj.k     = linspace(0, k_max_psd, Nk_psd);
            psd_vals  = zeros(1, Nk_psd);

            if nz == 1
                % 2D Hankel transform
                for ik = 1:Nk_psd
                    kr = obj.k(ik) * r_norm;
                    J0 = besselj(0, kr);
                    psd_vals(ik) = (1/(2*pi)) * trapz(r_norm, r_norm .* Racf .* J0);
                end
            else
                % 3D spherical Fourier-Bessel (sinc) transform
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
            end
            psd_vals = max(psd_vals, 0);

            obj.Phi = @(z) interp1(obj.k, psd_vals, z, 'makima', 0);
            out     = obj.Phi;

            r_plot  = r_norm(r_norm <= 6);
            S2_plot = S2_radial(r_norm <= 6);
            R_plot  = Racf(r_norm <= 6);
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
            
            % Use a logarithmic grid for integration to capture decay
            if(isempty(obj.r))
                r = logspace(-4, 3, 4096); 
            else
                r = obj.r;
            end
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
            r = logspace(-4, 3, 4096);  % real space grid
            R_vals = obj.R(r);
            
            % Wavenumber grid
            k = logspace(-3, 3, 4096);  % wavenumber grid
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
                        ix = mod(find(abs(xv - cx) <= R) - 1, nx) + 1;
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
                        ix = mod(find(abs(xv - cx) <= R) - 1, nx) + 1;
                        iy = mod(find(abs(yv - cy) <= R) - 1, ny) + 1;
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
                        ix = mod(find(abs(xv - cx) <= R) - 1, nx) + 1;
                        iy = mod(find(abs(yv - cy) <= R) - 1, ny) + 1;
                        iz = mod(find(abs(zv - cz) <= R) - 1, nz) + 1;
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

                % Diffusion coefficient m²/s (Eq. (5.12), Ryzhik et al, 1996)
                mat.Diffusivity = mat.v^2/(double(d)*(mat.Sigma-mat.Sigmapr));

                mat.meanFreeTime = 1/mat.Sigma;
                mat.meanFreePath = mat.v * mat.meanFreeTime;

                mat.transportMeanFreePath = double(d) * mat.Diffusivity / mat.v;
                mat.transportMeanFreeTime = mat.transportMeanFreePath/mat.v;

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
            if any(isnan(pdf))
                warning('There is a NaN inside the probability density function')
            end
            cdf = cumsum(pdf,"omitnan")*(pi/Nth);
            ind = find(diff(cdf)>1e-8);
            ind = unique([ind ind+1]);
            invcdf = griddedInterpolant(cdf(ind),xth(ind));
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
            %   DD       – object diameter (or length for 1-D)
            %   L        – box size: scalar (1-D), 1x2 (2-D), or 1x3 (3-D)

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
    end
end