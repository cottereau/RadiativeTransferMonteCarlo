function sigma = PSDF2sigma(mat,d)
% Compute the normalized differential scattering cross-sections based on
% statistics of the fluctuating material parameters of the wave equation
%
% d                            dimension of the problem
% material.acoustics           acoustics [true] or elastics [false]
% material.Frequency           readial frequency
% material.correlationLength   correlation length (assumed identical for
%                              all material parameters)
% material.spectralType        type of power spectrum density, assumed
%                              identical for all material parameters ('exp',
%                              'power_law', 'gaussian', 'triangular' or
%                              'low_pass'
% material.coefficients_of_variation 
%                              vector defining the coefficients of
%                              variation of random medium's characteristic
%                              parameters
%                              * acoustics: (1) compressibility (2) density
%                              * elastics: (1) lambda (2) mu (3) density
% material.correlation_coefficients
%                              vector defining the correlation coeffients
%                              of random medium's characteristic parameters
%                              * acoustics: (1) (compressibility,density)
%                              * elastics: (1) (lambda,mu), (2) (lambda,density),
%                                          (3) (mu,density)
% Output                       Differential scattering cross-sections expressed
%                              in units of inverse of time.

% The formulas are based on
% L. Rhyzik, G. Papanicolaou, J. B. Keller. Transport equations for elastic
% and other waves in random media. Wave Motion 24, pp. 327-370, 1996.
% doi: 10.1016/S0165-2125(96)00021-2

% In Ryzhik et al, the PSDFs are those of the fractional parts (normalized)
% of the corresponding random fields. The corr_matrix thus contains the
% coefficients of variation in diagonal and correlation coefficients in off-diagonal
% instead of the standard deviation and covariance, respectively.

% Following Khazaie et al 2017, PSDFs are supposed as factorizable. 
% Extending to the non-factorizable case is straightforward.

% constants
acoustics = mat.acoustics;
coeffs_variation = mat.coefficients_of_variation;
corrcoefs = mat.correlation_coefficients;
if any(abs(corrcoefs)>1)
    error('Absolute values of correlation coefficients should be less than 1!')
end

% The following normalized PSDF kernels are taken from Khazaie et al 2016
% warning: these are for 3D: formulas should depend on dimensionality
switch mat.spectralType
    case 'exp'
        S = @(z) 1./(8*pi^2*(1+z.^2/4).^2);
    case 'power_law'
        S = @(z) 1./(pi^4)*exp(-2*z/pi);
    case 'gaussian'
        S = @(z) 1./(8*pi^3)*exp(-z.^2/4/pi);
    case 'triangular'
        S = @(z) (3/8/pi^4)*(1-z/2/pi).*heaviside(2*pi-z);
    case 'low_pass'
        S = @(z) (2/9/pi^4)*heaviside(3*pi/2-z);
    otherwise
        disp(['Power spectrum density (mat.spectralType) should be' ...
            ' ''exp'',''power_law'',''gaussian'',''triangular'' or ''low_pass'''])
end

% acoustics 
if acoustics
    if d==2
        error(['Total scattering cross-sections for 2D case are to be ' ...
               'implemented in future releases of the code'])
    elseif d==3
        zeta = mat.Frequency/mat.v*mat.correlationLength;
        delta_kk = coeffs_variation(1); % coefficient of variation of compressibility
        delta_rr = coeffs_variation(2); % coefficient of variation of density
        rho_kr = corrcoefs; % correlation of compressibility and density

        % [Ryzhik et al, 1996; Eq. (1.3)]
        sigma = @(th) (pi/2)*zeta^3*(cos(th).^2*delta_rr^2 + ...
                       2*cos(th)*rho_kr*delta_kk*delta_rr + delta_kk^2) ...
                       .*S(zeta.*sqrt(2*(1-cos(th))))*mat.Frequency;
    else
        error('Dimension of the problem should be 2 or 3')
    end

    % elastics
else
    if d==2
        error(['Total scattering cross-sections for 2D case are to be ' ...
               'implemented in future releases of the code'])
    elseif d==3
        K = mat.vp/mat.vs;
        zetaP = mat.Frequency/mat.vp*mat.correlationLength;
        zetaS = K*zetaP;
        delta_ll = coeffs_variation(1); % squared coefficient of variation of lambda
        delta_mm = coeffs_variation(2); % squared coefficient of variation of mu
        delta_rr = coeffs_variation(3); % squared coefficient of variation of density
        rho_lm = corrcoefs(1); % correlation coefficient between lambda and mu
        rho_lr = corrcoefs(2); % correlation coefficient between lambda and density
        rho_mr = corrcoefs(3); % correlation coefficient between mu and density

        % [Ryzhik et al; Eq. (1.3)] and [Turner, 1998; Eq. (3)]
        sigmaPP = @(th) (pi/2)*zetaP^3* ...
            ( (1-2/K^2)^2*delta_ll^2 + 4*(1/K^2-2/K^4)*rho_lm*delta_ll*delta_mm*cos(th).^2 ...
            + (4/K^4)*delta_mm^2*cos(th).^4 + delta_rr^2*cos(th).^2 ...
            + 2*(1-2/K^2)*rho_lr*delta_ll*delta_rr*cos(th) ...
            + (4/K^2)*rho_mr*delta_mm*delta_rr*cos(th).^3 ) ...
            .*S(zetaP.*sqrt(2*(1-cos(th))))*mat.Frequency;

        %  Ryzhik et al; Eqs. (4.56), (1.20), (1.22)
        sigmaPS = @(th) (pi/2)*K*zetaP^3* ...
            ( K^2*delta_rr^2 + 4*delta_mm^2*cos(th).^2 + 4*K*rho_mr*delta_mm*delta_rr*cos(th) )...
              .*(1-cos(th).^2).*S(zetaP.*sqrt(1+K^2-2*K*cos(th)))*mat.Frequency;

        %  Ryzhik et al; Eqs. (4.56), (1.20), (1.21)
        sigmaSP = @(th) (pi/2/K^3)*zetaS^3* ...
            ( delta_rr^2 + (4/K^2)*delta_mm^2*cos(th).^2 + (4/K)*rho_mr*delta_mm*delta_rr*cos(th) ) ...
              .*(1-cos(th).^2).*S(zetaS.*sqrt(1+1/K^2-2/K*cos(th)))*mat.Frequency;

        % Ryzhik et al; Eq. (4.54)
        sigmaSS_TT = @(th) (pi/2)*zetaS^3*delta_rr^2*(1+cos(th).^2)...
                            .*S(zetaS.*sqrt(2*(1-cos(th))))*mat.Frequency;
        sigmaSS_GG = @(th) (pi/2)*zetaS^3*delta_mm^2*(4*cos(th).^4-3*cos(th).^2+1)...
                            .*S(zetaS.*sqrt(2*(1-cos(th))))*mat.Frequency;
        sigmaSS_GT = @(th) (pi/2)*zetaS^3*rho_mr*delta_mm*delta_rr*(2*cos(th).^3)...
                            .*S(zetaS.*sqrt(2*(1-cos(th))))*mat.Frequency;
    end

    sigmaSS = @(th) sigmaSS_TT(th) + sigmaSS_GG(th) + sigmaSS_GT(th);

    sigma = {sigmaPP,sigmaPS; ...
             sigmaSP,sigmaSS};
end
end

function h = heaviside(h)
h(h>0) = 1;
h(h==0) = 1/2;
h(h<0) = 0;
end