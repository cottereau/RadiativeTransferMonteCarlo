classdef Medium
%MEDIUM  Homogeneous isotropic elastic solid (beta > 0) or inviscid fluid (beta = 0).
%
%   m = lamb.Medium(rho, alpha, beta)      % solid
%   m = lamb.Medium(rho, alpha)            % fluid (beta = 0)
%
%   Properties : rho [kg/m^3], alpha [m/s] (P / sound speed), beta [m/s] (S speed)
%   Dependent  : mu, lambda (Lame constants), isFluid
%
%   Vertical wavenumber on the physical Riemann sheet (eq:branch):
%       nu = -1i*sqrt(km^2 - k^2)   ->   Re(nu) >= 0
%   valid both for real omega and for the complex frequency omega + 1i*eps (eps > 0).

    properties (SetAccess = immutable)
        rho   (1,1) double
        alpha (1,1) double
        beta  (1,1) double
    end
    properties (Dependent)
        mu
        lambda
        isFluid
    end

    methods
        function obj = Medium(rho, alpha, beta)
            arguments
                rho   (1,1) double {mustBePositive}
                alpha (1,1) double {mustBePositive}
                beta  (1,1) double {mustBeNonnegative} = 0
            end
            if beta >= alpha*sqrt(3)/2
                error('lamb:Medium:velocities', 'Need beta < (sqrt(3)/2)*alpha (Poisson ratio > -1).');
            end
            obj.rho = rho; obj.alpha = alpha; obj.beta = beta;
        end

        function v = get.mu(obj),      v = obj.rho*obj.beta^2;             end
        function v = get.lambda(obj),  v = obj.rho*obj.alpha^2 - 2*obj.mu; end
        function v = get.isFluid(obj), v = obj.beta == 0;                  end

        function km = wavenumber(obj, omega, wave)
            %WAVENUMBER  k_m = omega/c_m for wave 'P' or 'S' (also accepts 'SV', 'SH').
            if upper(wave(1)) == 'P'
                km = omega/obj.alpha;
            elseif obj.isFluid
                error('lamb:Medium:fluidS', 'A fluid has no S wave.');
            else
                km = omega/obj.beta;
            end
        end

        function nu = verticalWavenumber(obj, k, omega, wave)
            %VERTICALWAVENUMBER  nu_m = sqrt(k^2 - k_m^2) with Re(nu) >= 0 (element-wise).
            km = obj.wavenumber(omega, wave);
            nu = -1i*sqrt(km.^2 - k.^2);
        end
    end
end
