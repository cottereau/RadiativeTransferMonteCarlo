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
            %   newobj = calcSigmaIsotropic (  );
            %
            % Inputs:
            %
            % Output:

            % The formulas are based on
            % L. Rhyzik, G. Papanicolaou, J. B. Keller. Transport equations for elastic
            % and other waves in random media. Wave Motion 24, pp. 327-370, 1996.
            % doi: 10.1016/S0165-2125(96)00021-2
            %it works only for 3D cases
            switch obj.acoustics
                case 1
                    % acoustics
                    zeta = (2*pi*obj.Frequency)/obj.v*obj.CorrelationLength;
                    delta_kk = obj.coefficients_of_variation(1); % coefficient of variation of compressibility
                    delta_rr = obj.coefficients_of_variation(2); % coefficient of variation of density
                    rho_kr = obj.correlation_coefficients; % correlation of compressibility and density

                    % [Ryzhik et al, 1996; Eq. (1.3)]
                    obj.sigma = {@(th) (pi/2)*zeta^3*(cos(th).^2*delta_rr^2 + ...
                        2*cos(th)*rho_kr*delta_kk*delta_rr + delta_kk^2) ...
                        .*obj.Phi(zeta.*sqrt(2*(1-cos(th))))*(2*pi*obj.Frequency)};

                case 0
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
            switch obj.d
                case 1
                    error(['Normalized PSDF kernels for 2D case are to be ' ...
                        'implemented in future releases of the code'])
                case 2
                    error(['Normalized PSDF kernels for 2D case are to be ' ...
                        'implemented in future releases of the code'])
                case 3
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
        function out = MonoDisperseSphere(obj, eta, D )
            %% MonoDisperseDisk
            % Compute the normalized power spectral density function for
            % mono disperse disk
            %
            % Syntax:
            %   MonoDisperseDisk (  );
            %
            % Inputs:
            %   eta : volume fractions
            %   D   : sphere diameter
            %   r   : radius vector
            %
            % Output:

            % Torquato S. and Stell G. (1985). Microstructure of two-phase
            % random media. V. The n-point matrix probability functions for
            % impenetrable spheres, J. Chem. Phys., 82(2), pp. 980-987.

            dim = 3;

            % discretization in Fourier space
            obj.k = linspace(0,6/D,4094);
            % discretization in space
            raux = linspace(0,5*D,4096);

            % number density of spheres
            rhoS = 3*eta/(4*pi);
            % Verlet Weiss correction
            %eta = eta - (eta^2)/16;
            %D = (6*eta/pi/rho)^(1/3)

            % normalization of the radii
            raux = 2*raux/D;

            % Union volume of two spheres of unit radius
            switch dim
                case 1
                    V2 = 4*ones( size(raux) );
                    V2( raux<2 ) = 2+raux(raux<2);
                case 2
                    % V2 = 2*pi*ones( size(r) );
                    error('not implemented yet')
                case 3
                    V2 = 8*pi/3*ones( size(raux) );
                    V2( raux<2 ) = 4*pi/3*(1+3/4*raux(raux<2)-raux(raux<2).^3/16);
            end
            % Fourier transform of the direct correlation function
            % using the Percus-Yevick approximation
            if (dim==1)||(dim==2)
                error('not implemented yet')
            end
            l1 = (1+2*eta)^2/(1-eta)^4;
            l2 = -(1+eta/2)^2/(1-eta)^4;

            c = -4*pi./(obj.k.^3) .* ( l1*(sin(2*obj.k)-2*obj.k.*cos(2*obj.k)) + ...
                3*eta*l2./obj.k.*( 4*obj.k.*sin(2*obj.k) + (2-4*obj.k.^2).* ...
                cos(2*obj.k) - 2 ) + ...
                eta*l1./(2*obj.k.^3) .* (( 6*obj.k.^2 - 3 -2*obj.k.^4 ) .* ...
                cos(2*obj.k) + ...
                (4*obj.k.^3-6*obj.k).*sin(2*obj.k) + 3 ));
            c(1) = -8*pi/3*((4+eta)*l1+18*eta*l2);

            % % Verlet L and Weis J
            % lambda1 = (1+2*eta)/(1-eta)^4;
            % lambda2 = -(1+0.5*eta)^2/(1-eta)^4;
            % c11 = -lambda1-6*eta*lambda2.*(r/2)-0.5*eta*lambda1*(r/2).^2;
            % c11(r>1)=0;
            % c1_aux=fft(c11,2*size(k,2));
            % c1 = c1_aux(1:size(k,2));

            % Fourier transform of the total correlation function
            % using the Ornstein-Zernike relation
            h = c./(1-rhoS*c);
            % h = c1./(1-rho*c1);
            % Fourier transform of the Heaviside function
            switch dim
                case 1
                    m = 2 * sin(obj.k)./obj.k;
                    m(1) = 2;
                case 2
                    m = 2*pi * besselj(1,obj.k)./obj.k;
                    m(1) = pi;
                case 3
                    m = 4*pi * ( ((sin(obj.k)./obj.k) - cos(obj.k))./(obj.k.^2) );
                    m(1) = 4*pi/3;
            end
            % computation of 2-point matrix probability function
            M = zeros( size(raux) );
            if (dim==1)||(dim==2)
                error('not implemented yet')
            end
            M(1) = -(eta/rhoS)^2;
            for i1 = 2:length(raux)
                M(i1) = 1/2/pi^2/raux(i1)*trapz( obj.k, h.*m.^2.*obj.k.*sin(obj.k*raux(i1)) );
            end
            S2 = 1 - rhoS*V2 + rhoS^2*M + eta^2;
            %trapz(r,S2)
            % computation of the particle autocorrelation function
            Racf = (S2 -(1-eta)^2) / eta / (1-eta);
            raux = D*raux/2;
            right_window_length = 200;
            hann_win = hann(2 * right_window_length);
            Racf(end-right_window_length+1:end) = Racf(end-right_window_length+1:end) .* hann_win(right_window_length+1:end)';

            obj.R = @(z)interp1(raux,Racf,z,'makima',0);
            LL = 2*integral(obj.R,0,inf);
            obj.CorrelationLength = LL;
            raux = raux/LL;
            obj.R = @(z)interp1(raux,Racf,z,'makima',0);
            dk = 1/mean(diff(raux));
            %obj.k = linspace(0,dk,numel(raux));
            obj.k = obj.k*obj.CorrelationLength;
            phi = zeros(1,length(obj.k));
            for i1 = 1:length(obj.k)
                phi(i1) = 1/2/pi^2*trapz( raux, sinc(raux * obj.k(i1)) .* raux.^2 .* Racf );
            end

            P=abs(phi.*conj(phi));

            figure
            subplot(2,1,1)
            plot(raux,Racf)
            subplot(2,1,2)
            plot(obj.k,P)

            obj.Phi = @(z)interp1(obj.k,P,z,'makima',0);
            %obj.RalcLc;
            %obj.k = obj.k/obj.CorrelationLength;

            out = obj.Phi;
        end
        %% function to evaluate PSDF: an image
        function out = GetPSDFromImage(obj,Im,dim,dx,dy)
            %% MonoDispersGetPSDFromImageeDisk
            % Compute the normalized power spectral density function from
            % an image
            %
            % Syntax:
            %   GetPSDFromImage (  );
            %
            % Inputs:
            %
            % Output:


            % check dx and dy
            if ~exist('dx','var')
                dx = 1;
            end
            if ~exist('dy','var')
                dy = 1;
            end

            Freq = [1/dx 1/dy];
            x = 0:dx:size(Im,dim)*dx-dx;
            %
            %https://mathworld.wolfram.com/Wiener-KhinchinTheorem.html

            %convert to gray
            Im = im2gray(Im);
            %convert to double
            Im=double(Im)/256;

            %calc correlation length
            kappa = obj.v^2*1;
            variance = kappa*obj.coefficients_of_variation(1);

            [rcor,lags] = xcorr(mean(Im,dim),'normalized');
            x = lags*dx;
            %rcor = rcor/size(Im,dim);

            Lc = trapz(x,rcor)/variance;
            obj.CorrelationLength = Lc;
            obj.r = x/Lc;
            ind = obj.r<0;
            rcor(ind) = [];
            obj.r(ind) = [];
            obj.R = @(z)interp1(obj.r,rcor,z,'makima',0);

            dr = 1/mean(diff(obj.r));
            obj.k = linspace(0,dr/50,numel(rcor));

            phi = zeros(1,length(obj.k));
            for i1 = 1:length(obj.k)
                phi(i1) = 1/2/pi^2*trapz( obj.r, sinc(obj.r * obj.k(i1)) .* obj.r.^2 .* rcor );
            end
            obj.Phi = @(z)interp1(obj.k,phi,z,'makima',0);
            out = obj.Phi;
        end
        %% correlation length
        function Lc = CalcLc(obj)
            %2/variance integral( obj.R)| (R+ )
            % assuing rho = 1
            %kappa = obj.v^2*1;
            %variance = kappa*obj.Roefficients_of_variation(1);
            Lc = integral(obj.R,0,Inf)*2;
            %obj.CorrelationLength = Lc;
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
                plot(z,obj.sigma{1}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('Differential Scattering Cross-Sections [-]')
                grid on
                box on
                set(gca,'FontSize',14)
            else
                subplot(2,2,1)
                plot(z,obj.sigma{1,1}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('P2P [-]')
                grid on
                box on
                set(gca,'FontSize',14)

                subplot(2,2,2)
                plot(z,obj.sigma{1,2}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('P2S [-]')
                grid on
                box on
                set(gca,'FontSize',14) 

                subplot(2,2,3)
                plot(z,obj.sigma{2,1}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('S2P [-]')
                grid on
                box on
                set(gca,'FontSize',14)

                subplot(2,2,4)
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

                
                %xlabel('Angle [rad]')
                title('Differential Scattering Cross-Sections')
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
    end
end