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
        d               int8    = 3;% dimension
        %mat             struct  = struct.empty % PSDF parameter
        type            char    = 'isotropic'
        acoustics        

        v
        vp
        vs
        Frequency
        coefficients_of_variation
        correlation_coefficients

        SpectralLaw     char    = '';% PSDF name
        SpectralParam   struct  = struct.empty % PSDF parameter

        CorrelationLength       = [];% correlation length


        sigma           cell    = cell.empty; %Differential Scattering Cross-Sections

        Sigma
        Sigmapr
        invcdf
        D
        meanFreeTime
        P2P
        S2S

        Phi                     = []; %power spectral density / function_handle
        k               double  = [];% wavenumber vector
        R                       = []; %Correlation / function_handle
        r               double  = [];% r vector
        timeSteps               = 0;% time Steps : 0=small 1=large

    end
    properties (SetAccess = private, Hidden = true)
        Type_def = {'isotropic'}; %the anisotropic should be implemented
        SpectralLaw_def = {'exp','power_law','gaussian','triangular','low_pass','VonKarman','monodispersedisk','monodisperseelipse','polydispersedisk','polydisperseelipse','image'};
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
                obj.getSPDF;
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
        function getSPDF(obj)
            %% getSPDF
            % Compute the normalized power spectral density function to use
            % it in the differential scattering cross-section
            %
            % Syntax:
            %   getSPDF (  );
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
                case 2
                    error(['Normalized PSDF kernels for 2D case are to be ' ...
                        'implemented in future releases of the code'])
                case 3
                    switch obj.SpectralLaw
                        case 'exp'
                            %obj.SpectralParam = struct('Lc',obj.CorrelationLength);
                            obj.Exponential(obj.CorrelationLength);
                        case 'power_law'
                            %obj.SpectralParam = struct('Lc',obj.CorrelationLength);
                            obj.PowerLaw(obj.CorrelationLength);
                        case 'gaussian'
                            %obj.SpectralParam = struct('Lc',obj.CorrelationLength);
                            obj.Gaussian(obj.CorrelationLength);
                        case 'triangular'
                            %obj.SpectralParam = struct('Lc',obj.CorrelationLength);
                            obj.Triangular(obj.CorrelationLength);
                        case 'low_pass'
                            %obj.SpectralParam = struct('Lc',obj.CorrelationLength);
                            obj.LowPass(obj.CorrelationLength);
                        case 'VonKarman'
                            if ~isfield(obj.SpectralParam,'H') && ~isfield(obj.SpectralParam,'nu')
                                error('Please add the VonKarman parameters for the PSDF')
                            end
                            obj.VonKarman(obj.SpectralParam);
                        case 'monodispersedisk'
                            if ~isfield(obj.SpectralParam,'rho') && ~isfield(obj.SpectralParam,'d')
                                error('Please add the mono disperse disk parameters for the PSDF')
                            end
                            obj.MonoDisperseDisk(obj.SpectralParam);
                        case 'monodisperseelipse'
                            if ~isfield(obj.SpectralParam,'rho') && ~isfield(obj.SpectralParam,'a') && ~isfield(obj.SpectralParam,'b')
                                error('Please add the mono disperse disk parameters for the PSDF')
                            end
                            obj.MonoDisperseElipse(obj.SpectralParam);
                        case 'polydispersedisk'
                            if ~isfield(obj.SpectralParam,'rho')
                                error('not defined yet')
                            end
                            obj.PolyDisperseDisk(obj.SpectralParam);
                        case 'polydisperseelipse'
                            if ~isfield(obj.SpectralParam,'rho')
                                error('not defined yet')
                            end
                            obj.PolyDisperseElipse(obj.SpectralParam);
                        case 'image'
                            if ~isfield(obj.SpectralParam,'ImagePath') && ~isfield(obj.SpectralParam,'dx') &&~isfield(obj.SpectralParam,'dy')
                                error('not defined yet')
                            end
                            obj.GetPSDFromImage(obj.SpectralParam);
                    end
            end
        end
        %% function to evaluate PSDF:
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

            out = @(z) 1./(8*pi^2*(1+z.^2/4).^2);
            obj.R = @(z) exp(-2*z);
            obj.Phi = out;
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
            out = @(z) 1./(pi^4)*exp(-2*z/pi);
            obj.Phi = out;
            obj.R = @(z) 1./(1+(pi^2*z.^2/4))^2;
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
            out = @(z) 1./(8*pi^3)*exp(-z.^2/4/pi);
            obj.Phi = out;
            obj.R = @(z)exp(-pi*z.^2);
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
        function out = VonKarman(obj,nu,H)
            %% VonKarman
            % Compute the normalized power spectral density function for
            % Von Karman
            %
            % Syntax:
            %   VonKarman (  );
            %
            % Inputs:
            %   nu : Hurst number
            %   H  : ?
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
            obj.SpectralParam = struct('nu',nu,'H',H);
            obj.SpectralLaw = 'VonKarman';

            error('Not fully implemented')
            r = sqrt(x^2/ax^2+y^2/ay^2+z^2/az^2);
            k=sqrt(kx^2*ax^2+ky^2*ay^2+kz^2*az^2);
            bessel = besseli(nu,r);
            C = 4*pi*nu*H^2*r^nu*bessel/bessel(1);
            S = (4*pi*nu*H^2/bessel(1))*(ax^2+ay^2+az^2)/((1+k^2^(nu+1.5)));
            out = @(z) 1;
            obj.Phi = out;
            obj.R = out;

            obj.CorrelationLength = obj.RalcLc;

        end
        function out = MonoDisperseDisk(obj, eta, D )
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
            obj.k = linspace(0,4/D,4094);
            % discretization in space
            r = linspace(0,10*D,4096);

            % number density of spheres
            rho = 3*eta/(4*pi);
            % Verlet Weiss correction
            %eta = eta - (eta^2)/16;
            %D = (6*eta/pi/rho)^(1/3)

            % normalization of the radii
            r = 2*r/D;

            % Union volume of two spheres of unit radius
            switch dim
                case 1
                    V2 = 4*ones( size(r) );
                    V2( r<2 ) = 2+r(r<2);
                case 2
                    % V2 = 2*pi*ones( size(r) );
                    error('not implemented yet')
                case 3
                    V2 = 8*pi/3*ones( size(r) );
                    V2( r<2 ) = 4*pi/3*(1+3/4*r(r<2)-r(r<2).^3/16);
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
            h = c./(1-rho*c);
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
            M = zeros( size(r) );
            if (dim==1)||(dim==2)
                error('not implemented yet')
            end
            M(1) = -(eta/rho)^2;
            for i1 = 2:length(r)
                M(i1) = 1/2/pi^2/r(i1)*trapz( obj.k, h.*m.^2.*obj.k.*sin(obj.k*r(i1)) );
            end
            S2 = 1 - rho*V2 + rho^2*M + eta^2;
            trapz(r,S2)
            % computation of the particle autocorrelation function
            R = (S2 -(1-eta)^2) / eta / (1-eta);

            obj.R = @(z)interp1(r,R,z,'makima',0);
            LL = 2*integral(obj.R,0,inf);
            obj.CorrelationLength = LL;
            %obj.k = obj.k*LL;
            phi = zeros(1,length(obj.k));
            for i1 = 1:length(obj.k)
                phi(i1) = 1/2/pi^2*trapz( r, sinc(r * obj.k(i1)) .* r.^2 .* R );
            end
           
            P=abs(((phi.*conj(phi))));

            figure
            subplot(2,1,1)
            plot(r/LL,R)
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
        function PlotPSD(obj)
            figure
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
        function PlotCorrelation(obj)
            figure
            x = linspace(0,10,2048);
            plot(x,obj.R(x),'LineWidth',2)
            xlabel('Normalized lag distance [-]')
            ylabel('Correlation [-]')
            grid on
            box on
            set(gca,'FontSize',14)
        end
        function plotsigma(obj)
            figure
            z = linspace(0,2*pi,2*2048);
            if obj.acoustics
                plot(z,obj.sigma{1}(z),'LineWidth',2)
                xlabel('Angle [rad]')
                ylabel('Differential Scattering Cross-Sections [-]')
                grid on
                box on
                set(gca,'FontSize',14)
            else
                error('Not implemented')
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
                error('Not implemented')
            end
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
    end
end