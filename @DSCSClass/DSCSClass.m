%% DSCSClass
% Class to deal with Differential Scattering Cross-Sections (DSCS)
%
%
% DSCSClass ()
%
% See also
%
%
classdef DSCSClass < handle
    properties
        % Required Properties
        d               int8
        mat             struct  = struct.empty % PSDF parameter
        type            char    = 'isotropic'

        SpectralLaw     char    = '';% PSDF name
        SpectralParam   struct  = struct.empty % PSDF parameter

        sigma           cell    = cell.empty; %Differential Scattering Cross-Sections
        
        S                       = []; %power spectral density / function_handle
        C                       = []; %Correlation / function_handle

    end
    properties (SetAccess = private, Hidden = true)
        Type_def = {'isotropic'}; %the anisotropic should be implemented
        SpectralLaw_def = {'exp','power_law','gaussian','triangular','low_pass','VonKarman','monodispersedisk','monodisperseelipse','polydispersedisk','polydisperseelipse','image'};
    end
    methods
        function obj = DSCSClass(mat,d)
            %% DSCSClass
            % DSCSClass contructor
            %
            % Syntax:
            %   newobj = DSCSClass (  );
            %
            % Inputs:
            %     mat : material structure
            %     d   : dimension of the problem

            if d~=3
                warning(['Total scattering cross-sections for 2D case are to be ' ...
                    'implemented in future releases of the code'])
            end

            obj.mat = mat;
            obj.d = d;
            obj.SpectralLaw = mat.spectralType;
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
            if numel(obj) == 1
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
            if isempty(obj.S)
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
            switch obj.mat.acoustics
                case 1
                    % acoustics
                    zeta = obj.mat.Frequency/obj.mat.v*obj.mat.correlationLength;
                    delta_kk = obj.mat.coefficients_of_variation(1); % coefficient of variation of compressibility
                    delta_rr = obj.mat.coefficients_of_variation(2); % coefficient of variation of density
                    rho_kr = obj.mat.correlation_coefficients; % correlation of compressibility and density

                    % [Ryzhik et al, 1996; Eq. (1.3)]
                    obj.sigma = {@(th) (pi/2)*zeta^3*(cos(th).^2*delta_rr^2 + ...
                        2*cos(th)*rho_kr*delta_kk*delta_rr + delta_kk^2) ...
                        .*obj.S(zeta.*sqrt(2*(1-cos(th))))*obj.mat.Frequency};
                case 0
                    %elastic
                    K = obj.mat.vp/obj.mat.vs;
                    zetaP = obj.mat.Frequency/obj.mat.vp*obj.mat.correlationLength;
                    zetaS = K*zetaP;
                    delta_ll = obj.mat.coefficients_of_variation(1); % squared coefficient of variation of lambda
                    delta_mm = obj.mat.coefficients_of_variation(2); % squared coefficient of variation of mu
                    delta_rr = obj.mat.coefficients_of_variation(3); % squared coefficient of variation of density
                    rho_lm = obj.mat.correlation_coefficients(1); % correlation coefficient between lambda and mu
                    rho_lr = obj.mat.correlation_coefficients(2); % correlation coefficient between lambda and density
                    rho_mr = obj.mat.correlation_coefficients(3); % correlation coefficient between mu and density

                    % [Ryzhik et al; Eq. (1.3)] and [Turner, 1998; Eq. (3)]
                    sigmaPP = @(th) (pi/2)*zetaP^3* ...
                        ( (1-2/K^2)^2*delta_ll^2 + 4*(1/K^2-2/K^4)*rho_lm*delta_ll*delta_mm*cos(th).^2 ...
                        + (4/K^4)*delta_mm^2*cos(th).^4 + delta_rr^2*cos(th).^2 ...
                        + 2*(1-2/K^2)*rho_lr*delta_ll*delta_rr*cos(th) ...
                        + (4/K^2)*rho_mr*delta_mm*delta_rr*cos(th).^3 ) ...
                        .*obj.S(zetaP.*sqrt(2*(1-cos(th))))*obj.mat.Frequency;

                    %  Ryzhik et al; Eqs. (4.56), (1.20), (1.22)
                    sigmaPS = @(th) (pi/2)*K*zetaP^3* ...
                        ( K^2*delta_rr^2 + 4*delta_mm^2*cos(th).^2 + 4*K*rho_mr*delta_mm*delta_rr*cos(th) )...
                        .*(1-cos(th).^2).*obj.S(zetaP.*sqrt(1+K^2-2*K*cos(th)))*obj.mat.Frequency;

                    %  Ryzhik et al; Eqs. (4.56), (1.20), (1.21)
                    sigmaSP = @(th) (pi/4/K^3)*zetaS^3* ...
                        ( delta_rr^2 + (4/K^2)*delta_mm^2*cos(th).^2 + (4/K)*rho_mr*delta_mm*delta_rr*cos(th) ) ...
                        .*(1-cos(th).^2).*obj.S(zetaS.*sqrt(1+1/K^2-2/K*cos(th)))*obj.mat.Frequency;

                    % Ryzhik et al; Eq. (4.54)
                    sigmaSS_TT = @(th) (pi/4)*zetaS^3*delta_rr^2*(1+cos(th).^2)...
                        .*obj.S(zetaS.*sqrt(2*(1-cos(th))))*obj.mat.Frequency;
                    sigmaSS_GG = @(th) (pi/4)*zetaS^3*delta_mm^2*(4*cos(th).^4-3*cos(th).^2+1)...
                        .*obj.S(zetaS.*sqrt(2*(1-cos(th))))*obj.mat.Frequency;
                    sigmaSS_GT = @(th) (pi/4)*zetaS^3*rho_mr*delta_mm*delta_rr*(4*cos(th).^3)...
                        .*obj.S(zetaS.*sqrt(2*(1-cos(th))))*obj.mat.Frequency;


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

            if any(abs(obj.mat.correlation_coefficients)>1)
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
                            if ~isfield(obj.mat,'correlationLength')
                                error('Please add the exponential parameters in material strucutre - correlationLength')
                            end
                            obj.SpectralParam = struct('Lc',obj.mat.correlationLength);
                            obj.Exponential;
                        case 'power_law'
                            if ~isfield(obj.mat,'correlationLength')
                                error('Please add the power law parameters in material strucutre - correlationLength')
                            end
                            obj.SpectralParam = struct('Lc',obj.mat.correlationLength);
                            obj.PowerLaw;
                        case 'gaussian'
                            if ~isfield(obj.mat,'correlationLength')
                                error('Please add the gaussian parameters in material strucutre - correlationLength')
                            end
                            obj.SpectralParam = struct('Lc',obj.mat.correlationLength);
                            obj.Gaussian;
                        case 'triangular'
                            if ~isfield(obj.mat,'correlationLength')
                                error('Please add the triangular parameters in material strucutre - correlationLength')
                            end
                            obj.SpectralParam = struct('Lc',obj.mat.correlationLength);
                            obj.Triangular;
                        case 'low_pass'
                            if ~isfield(obj.mat,'correlationLength')
                                error('Please add the low pass parameters in material strucutre - correlationLength')
                            end
                            obj.SpectralParam = struct('Lc',obj.mat.correlationLength);
                            obj.LowPass;
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
        function out = Exponential(obj)
            %% Exponential
            % Compute the normalized power spectral density function for
            % Exponential
            %
            % Syntax:
            %   Exponential (  );
            %
            % Inputs:
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            out = @(z) 1./(8*pi^2*(1+z.^2/4).^2);
            obj.S = out;
            obj.C = out;
        end
        function out = PowerLaw(obj)
            %% PowerLaw
            % Compute the normalized power spectral density function for
            % Power Law
            %
            % Syntax:
            %   PowerLaw (  );
            %
            % Inputs:
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            out = @(z) 1./(pi^4)*exp(-2*z/pi);
            obj.S = out;
            obj.C = out;
        end
        function out = Gaussian(obj)
            %% Gaussian
            % Compute the normalized power spectral density function for
            % Gaussian
            %
            % Syntax:
            %   Gaussian (  );
            %
            % Inputs:
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            out = @(z) 1./(8*pi^3)*exp(-z.^2/4/pi);
            obj.S = out;
            obj.C = out;
        end
        function out = Triangular(obj)
            %% Triangular
            % Compute the normalized power spectral density function for
            % Triangular
            %
            % Syntax:
            %   Triangular (  );
            %
            % Inputs:
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties
            out = @(z) (3/8/pi^4)*(1-z/2/pi).*obj.heaviside(2*pi-z);
            obj.S = out;
            obj.C = out;
        end
        function out = LowPass(obj)
            %% LowPass
            % Compute the normalized power spectral density function for
            % LowPass
            %
            % Syntax:
            %   LowPass (  );
            %
            % Inputs:
            %
            % Output:
            % The following normalized PSDF kernels are taken from
            % Khazaie et al 2016 - Influence of the spatial correlation
            % structure of an elastic random medium on its
            % scattering properties

            out = @(z) (2/9/pi^4)*obj.heaviside(3*pi/2-z);
            obj.S = out;
            obj.C = out;
        end
        function out = VonKarman(obj)
            %% VonKarman
            % Compute the normalized power spectral density function for
            % Von Karman
            %
            % Syntax:
            %   VonKarman (  );
            %
            % Inputs:
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
            
            r = sqrt(x^2/ax^2+y^2/ay^2+z^2/az^2);
            k=sqrt(kx^2*ax^2+ky^2*ay^2+kz^2*az^2);
            bessel = besseli(nu,r);
            C = 4*pi*nu*H^2*r^nu*bessel/bessel(1);
            S = (4*pi*nu*H^2/bessel(1))*(ax^2+ay^2+az^2)/((1+k^2^(nu+1.5)));
            out = @(z) 1;
            obj.S = out;
            obj.C = out;
        end
        function out = MonoDisperseDisk(obj)
            %% MonoDisperseDisk
            % Compute the normalized power spectral density function for
            % mono disperse disk
            %
            % Syntax:
            %   MonoDisperseDisk (  );
            %
            % Inputs:
            %
            % Output:

            % Torquato S. and Stell G. (1985). Microstructure of two-phase
            % random media. V. The n-point matrix probability functions for
            % impenetrable spheres, J. Chem. Phys., 82(2), pp. 980-987.


            % constants
            rho = 3*phi/(4*pi);
            % k = 0:0.01:100;

            pad_x = numel(r);
            x_fa = 1/mean(diff(r));
            k = 0:x_fa/pad_x:x_fa-x_fa/pad_x;

            eta = 4*pi*rho/3;
            l1 = (1+2*eta)^2/(1-eta)^4;
            l2 = -(1+eta/2)^2/(1-eta)^4;
            V2 = 8*pi/3*ones( size(r) );
            V2( r<2 ) = 4*pi/3*(1+3/4*r(r<2)-r(r<2).^3/16);
            % computation
            c = -4*pi./(8*k.^3) .* ( l1*(sin(2*k)-2*k.*cos(2*k)) + ...
                3*eta*l2./k.*( 4*k.*sin(2*k) + (2-4*k.^2).*cos(2*k) - 2 ) + ...
                eta*l1./(16*k.^3) .* (( 48*k.^2 - 24 -16*k.^4 ) .* cos(2*k) + ...
                (32*k.^3-48*k).*sin(2*k) + 24 ));
            c(1) = -pi/3*((4+eta)*l1+18*eta*l2);

            % c = -4*pi./(k.^3) .* ( l1*(sin(2*k)-2*k.*cos(2*k)) + ...
            % 3*eta*l2./k.*( 4*k.*sin(2*k) + (2-4*k.^2).* ...
            % cos(2*k) - 2 ) + ...
            % eta*l1./(2*k.^3) .* (( 6*k.^2 - 3 -2*k.^4 ) .* ...
            % cos(2*k) + ...
            % (4*k.^3-6*k).*sin(2*k) + 3 ));
            % c(1) = -8*pi/3*((4+eta)*l1+18*eta*l2);

            m = 4*pi./k .*( sin(k)./(k.^2) - cos(k)./k );
            m(1) = 4*pi/3;
            M = zeros( size(r) );
            for i1 = 2:length(r)
                M(i1) = 1/2/pi^2/r(i1) * ...
                    trapz( k, c./(1-rho*c).*m.^2.*k.*sin(k*r(i1)) )+ 16*pi^2/9;
            end
            C = (1 - rho*V2 + rho^2*M -(1-phi)^2 ) / phi / (1-phi);
            % C = (1 - rho*V2 + rho^2*M + eta^2 -(1-eta)^2 ) / eta / (1-eta);
            % C = (1 - rho*V2 + rho^2*M -(1-phi)^2 ) / phi / (1-phi);

            %
            %
            S = 1 - rho*V2 + rho^2*M ;
            out = @(z) 1;
            obj.S = out;
            obj.C = out;
        end
        function out = MonoDisperseElipse(obj)
            %% MonoDisperseElipse
            % Compute the normalized power spectral density function for
            % mono disperse elipsoid
            %
            % Syntax:
            %   MonoDisperseElipse (  );
            %
            % Inputs:
            %
            % Output:

            % Torquato S. and Stell G. (1985). Microstructure of two-phase
            % random media. V. The n-point matrix probability functions for
            % impenetrable spheres, J. Chem. Phys., 82(2), pp. 980-987.

            % constants
            rho = 3*phi/(4*pi);
            % k = 0:0.01:100;

            pad_x = numel(r);
            x_fa = 1/mean(diff(r));
            k = 0:x_fa/pad_x:x_fa-x_fa/pad_x;

            eta = 4*pi*rho/3;
            l1 = (1+2*eta)^2/(1-eta)^4;
            l2 = -(1+eta/2)^2/(1-eta)^4;
            V2 = 8*pi/3*ones( size(r) );
            V2( r<2 ) = 4*pi/3*(1+3/4*r(r<2)-r(r<2).^3/16);
            % computation
            c = -4*pi./(8*k.^3) .* ( l1*(sin(2*k)-2*k.*cos(2*k)) + ...
                3*eta*l2./k.*( 4*k.*sin(2*k) + (2-4*k.^2).*cos(2*k) - 2 ) + ...
                eta*l1./(16*k.^3) .* (( 48*k.^2 - 24 -16*k.^4 ) .* cos(2*k) + ...
                (32*k.^3-48*k).*sin(2*k) + 24 ));
            c(1) = -pi/3*((4+eta)*l1+18*eta*l2);

            % c = -4*pi./(k.^3) .* ( l1*(sin(2*k)-2*k.*cos(2*k)) + ...
            % 3*eta*l2./k.*( 4*k.*sin(2*k) + (2-4*k.^2).* ...
            % cos(2*k) - 2 ) + ...
            % eta*l1./(2*k.^3) .* (( 6*k.^2 - 3 -2*k.^4 ) .* ...
            % cos(2*k) + ...
            % (4*k.^3-6*k).*sin(2*k) + 3 ));
            % c(1) = -8*pi/3*((4+eta)*l1+18*eta*l2);

            m = 4*pi./k .*( sin(k)./(k.^2) - cos(k)./k );
            m(1) = 4*pi/3;
            M = zeros( size(r) );
            for i1 = 2:length(r)
                M(i1) = 1/2/pi^2/r(i1) * ...
                    trapz( k, c./(1-rho*c).*m.^2.*k.*sin(k*r(i1)) )+ 16*pi^2/9;
            end
            C = (1 - rho*V2 + rho^2*M -(1-phi)^2 ) / phi / (1-phi);
            % C = (1 - rho*V2 + rho^2*M + eta^2 -(1-eta)^2 ) / eta / (1-eta);
            % C = (1 - rho*V2 + rho^2*M -(1-phi)^2 ) / phi / (1-phi);

            %
            %
            S = 1 - rho*V2 + rho^2*M ;

            out = @(z) 1;
            obj.S = out;
            obj.C = out;
        end
        function out = PolyDisperseDisk(obj)
            %% PolyDisperseDisk
            % Compute the normalized power spectral density function for
            % poly disperse disk
            %
            % Syntax:
            %   PolyDisperseDisk (  );
            %
            % Inputs:
            %
            % Output:

            out = @(z) 1;
            obj.S = out;
            obj.C = out;
        end
        function out = PolyDisperseElipse(obj)
            %% PolyDisperseElipse
            % Compute the normalized power spectral density function for
            % poly disperse elipsoid
            %
            % Syntax:
            %   PolyDisperseElipse (  );
            %
            % Inputs:
            %
            % Output:

            out = @(z) 1;
            obj.S = out;
            obj.C = out;
        end
        %% function to evaluate PSDF: an image
        function out = GetPSDFromImage(obj)
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

            Im=imread('rice.png');
            % check dx and dy
            if ~exist('dx','var')
                dx = 1;
            end
            if ~exist('dy','var')
                dy = 1;
            end

            %x = 0:dx:size(Im,1)-1;
            %y = 0:dy:size(Im,2)-1;
            kx = linspace(-1/dx/2,1/dx/2,size(Im,1));
            ky = linspace(-1/dy/2,1/dy/2,size(Im,2));

            %
            %https://mathworld.wolfram.com/Wiener-KhinchinTheorem.html

            %convert to double
            Im=double(Im);

            % 2D window
            r = 0.1;
            w1 = repmat(tukeywin(size(Im,1),r),1,size(Im,2));
            w2 = repmat(tukeywin(size(Im,2),r)',size(Im,1),1);
            Im = Im.*w1.*w2;

            %subtract mean
            Im=Im-mean(Im(:));

            %normalize magnitude
            Im=Im/sqrt(sum(Im(:).^2));

            %compute fft2
            I=fft2(Im);

            % auto correlation
            % A=real(fftshift(ifft2(I.*conj(I))));

            %Power
            P=abs(fftshift((I.*conj(I))))/numel(Im);
            %% calc the radius and theta
            rho = zeros(numel(kx),numel(ky));
            thetaI = zeros(numel(kx),numel(ky));
            for ix= 1 : numel(kx)
                for iy= 1 : numel(ky)
                    rho(ix,iy) = sqrt(kx(ix)^2+ky(iy)^2);
                    thetaI(ix,iy) = atan2(ky(iy),kx(ix));
                end
            end
            %% make the interpolation
            normk = linspace(0,max(rho(:)),min([numel(kx),numel(ky)]));
            out = zeros(size(normk));
            for i = 1:numel(normk)
                %out(i) = griddata(rho,thetaI,poissrnd,normk(i),theta,'cubic');
                out(i) = griddata(rho,thetaI,P,normk(i),theta,'nearest'); %#ok<GRIDD>
            end
            % return

            %Pinter = scatteredInterpolant(rho,thetaI,PI);
            %out = Pinter(normk,theta*ones(1,numel(normk)));

            %rho = reshape(rho,1,numel(rho))';
            %thetaI = reshape(rho,1,numel(thetaI))';
            %PI = reshape(P,1,numel(P))';


            out = @(z) 1;
            obj.S = out;
            obj.C = out;
        end
        %% correlation length
        function Lc = CalcLc(obj)
            %2/variancia integral em R+ do obj.S
            % assuing rho = 1
            kappa = obj.mat.v^2*1;
            variance = kappa*obj.mat.coefficients_of_variation(1);
            Lc = integral(obj.S,0,Inf);
            Lc = Lc*2/variance;
        end
        %% plot
        function Plotsigma(obj)
            figure
            z = linspace(0,10,2*2048);
            plot(z,obj.S(z));
            xlabel('Normalized wavenumer [-]')
            ylabel('Power Spectral Density [-]')
            grid on
            box on
        end
        %% OTHER
        function h = heaviside(h)
            h(h>0) = 1;
            h(h==0) = 1/2;
            h(h<0) = 0;
        end
    end
end