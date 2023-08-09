function sigma = PSDF2sigma( d, material )
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
% material.correlationMatrix   matrix giving the variance and correlation
%                              between material parameters

% Note 1: sigma has the unit 1/[T], as such, here the diff. scat. cross-sections
% are normalized by the angular frequency w.

% Note 2: in Ryzhik, the PSDFs are those of the fractional parts (normalized)
% of the corresponding random fields. The corr_matrix thus contains the
% coefficients of variation in diagonal and correlation coefficients in off-diagonal
% instead of the standard deviation and covariance, respectively.

% Note 3: following Khazaie et al 2017, PSDFs are supposed as factorizable
% but extending to the non-factorizable case is straightforward

% Possible bugs to be checked : 
% 1) How sigmaSP and sigmaSS, which are 2*1 & 2*2
% matrices, change in 2D? In this version of the code, in 2D phi=0 and the
% operators are no longer proportional to the identity matrix (why?).
% 2) See how the multiplication by the identity or other matrices should be
% done in the code (it is not done everywhere).

% constants
acoustics = material.acoustics;
zeta = material.Frequency/material.v*material.correlationLength;
C = material.correlationMatrix;

% test the symmetry of the correlation matrix
if ~issymmetric(material.correlationMatrix)
    disp('correlationMatrix should be a symmetric matrix !')
end

% The following normalized PSDF kernels are taken from Khazaie et al 2016
% warning: these are for 3D: formulas should depend on dimensionality
switch material.spectralType
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
        disp(['Power spectrum density (material.spectralType) should be' ...
     ' ''exp'',''power_law'',''gaussian'',''triangular'' or ''low_pass'''])
end

% Excep for sigmaSS, the constants in 2D & 3D cross sections are respectively 
% pi/2 & pi^2. We can generalize in d-dimension using the following formula
coeff = pi^(3*d/2-2)/2/gamma(d/2);

% acoustics [Ryzhik et al, 1996 Eq. (1.3)]
if acoustics 
    std_kk = C(1,1); % variance of the compressibility
    std_rr = C(2,2); % variance of the density
    std_kr = C(1,2); % correlation of compressibility and density
    sigma = @(th) coeff*zeta^d*(cos(th).^2*std_rr^2 + 2*cos(th)*std_kr + std_kk^2) ...
        .*S(zeta.*sqrt(2*(1-cos(th)))).*sin(th).^(d-2);

% elastics
else
    K = material.vp/material.vs;
    zetaP = zeta;
    zetaS = K*zetaP;
    std_ll = C(1,1);
    std_mm = C(2,2);
    std_rr = C(3,3);
    std_lm = C(1,2);
    std_lr = C(1,3);
    std_mr = C(2,3);
    
    % [Ryzhik et al, 1996; Eq. (1.3)] and [Turner, 1998; Eq. (3)  
    sigmaPP = coeff*zetaP^d * ...
           ( (1-2/K^2)^2*std_ll^2 + 4*(1/K^2-2/K^4)*std_lm*cos(th).^2 + ...
             (4/K^4)*std_mm^2*cos(th).^2 + std_rr^2*cos(th).^2 + ...
             2*(1-2/K^2)*std_lr*cos(th) + (4/K^2)*std_mr*cos(th).^3 ) ...
             .*sin(th).^(d-2).*S(zetaP.*sqrt(2*(1-cos(th))));

    %  [Ryzhik, 1996; Eq. (4.56)] and [Khazaie, PhD thesis; Eq. (1.81)]
    sigmaPS = coeff*K*zetaP^d * ...
          ( K^2*std_rr^2 + 4*std_mm^2*cos(th).^2 + 4*K*std_mr*cos(th) )...
          .*sin(th).^(d-2).*S(zetaP.*sqrt(1+K^2-2*K*cos(th)));

    % See Ryzhik 1996 Eq. (1.19) & Khazaie 2015 (PhD thesis) Eq. (1.82)
    % In 3D, for each value of theta (th), sigmaSP should be multiplied by
    % the identity matrix, i.e. eye(2,2)
    % In 2D, it should be multiplied by [1  0 ; 0  0];
    sigmaSP = coeff/(2*K^4)*zetaS^d * ...
          ( std_rr^2 + (4/K^2)*std_mm^2*cos(th).^2 + (4/K)*std_mr*cos(th) ) ...
          .*(1-cos(th).^2).*sin(th).^(d-2).*S(zetaS.*sqrt(1+1/K^2-2/K*cos(th)));

    % Ryzhik Eq. (4.54) & Khazaie 2015 (PhD thesis) Eqs. (1.84) to (1.86) 
    % In 3D, for each value of theta (th), sigmaSS should be multiplied by
    % the identity matrix, i.e. eye(2,2)
    if d==3
        sigmaSS_TT = @(th) pi^(d-1)/2 * zetaS^d * std_rr^2*(1+cos(th).^2)...
            .*sin(th).^(d-2).*S(zetaS.*sqrt(2*(1-cos(th)))); 
        sigmaSS_GG = @(th) pi^(d-1)/2 * zetaS^d * std_mm^2*(4*cos(th).^4-3*cos(th).^2+1)...
            .*sin(th).^(d-2).*S(zetaS.*sqrt(2*(1-cos(th))));
        sigmaSS_GT = @(th) pi^(d-1)/2 * zetaS^d * std_mr*(2*cos(th).^3)...
            .*sin(th).^(d-2).*S(zetaS.*sqrt(2*(1-cos(th))));
    elseif d==2
        sigmaSS_TT = @(th) pi^(d-1)/2 * zetaS^d * std_rr^2 ...
            *[ cos(th).^2 zeros(size(th)); zeros(size(th)) ones(size(th))].*S(zetaS.*sqrt(2*(1-cos(th))));
        sigmaSS_GG = @(th) pi^(d-1)/2 * zetaS^d * std_mm^2 ...
            *[ cos(th).^2.*(2*cos(th).^2-1).^2 zeros(size(th)); zeros(size(th)) cos(th).^2].*S(zetaS.*sqrt(2*(1-cos(th))));
        sigmaSS_GT = @(th) pi^(d-1)/2 * zetaS^d * std_mr ...
            *[ cos(th).*(2*cos(th).^2-1).^2 zeros(size(th)); zeros(size(th)) cos(th)].*S(zetaS.*sqrt(2*(1-cos(th))));
    end

    sigmaSS = @(th) sigmaSS_TT(th) + sigmaSS_GG(th) + sigmaSS_GT(th);

    sigma = cell(sigmaPP,sigmaPS,sigmaSP,sigmaSS);

end
end

function h = heaviside(h)
h(h>0) = 1;
h(h==0) = 1/2;
h(h<0) = 0;
end