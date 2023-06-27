function [Sigma,invcdf] = prepareSigmaOne_beta(source, material, geometry)
freq = source.frequency; % This field should be added to the source struct in Main
v = material.v;
lc = material.lc; % Correlation distance => field to be added either to material or another separate variable in main
sigma = material.sigma; % Maybe separating the statistical parameters in a separate variable?
% Note : psdf (type of the PSDF) should be added as a field to either material or another variable
% Note : the covariance matrix, containing the correlations between different pairs of elastic properties should be further added as a field?
% Note : in material (or a new variable) define if the density is random or not (material.density_random)
d = geometry.dimension;

Nth = 1000;
xth = linspace(0,pi,Nth);
% The following normalized PSDF kernels are taken from Khazaie et al 2016 
if strcmpi(material.psdf,'exp')
    S = @(z) 1./(8*pi^2*(1+z.^2/4).^2);
elseif strcmpi(material.psdf,'power_law')
    S = @(z) 1./(pi^4)*exp(-2*z/pi);
elseif strcmpi(material.psdf,'gaussian')
    S = @(z) 1./(8*pi^3)*exp(-z.^2/4/pi);
elseif strcmpi(material.psdf,'triangular')
    S = @(z) (3/8/pi^4)*(1-z/2/pi).*heaviside(2*pi-z);
elseif strcmpi(material.psdf,'low_pass')
    S = @(z) (2/9/pi^4)*heaviside(3*pi/2-z);
else
    disp("PSDF should be 'exp','power_law','gaussian','triangular' or 'low_pass'")
end

% Normalized frequency
zeta = 2*pi*freq*lc/v;

% adimensional diff scattering cross-section ( := cross-section / v*L^(d-1) )
% Based on Ryzhik 1996
if material.acoustics && material.density_random
    % random density
    sigma = @(th,z) (pi/2)*z.*(cos(th)+1).^2.*S(z.*sqrt(2*(1-cos(th))));
elseif material.acoustics && ~material.density_random
    % deterministic density
    sigma = @(th,z) (pi/2)*z.*S(z.*sqrt(2*(1-cos(th))));
end

% adimensional total scattering cross-section ( := total scat cross_sec * L / v )
% Based on Ryzhik 1996
if material.acoustics && d==2
    Sigma_hat = @(z) integral(@(th) z.^2.*sigma(th,z),0,pi);
elseif material.acoustics && d==3
    Sigma_hat = @(z) integral(@(th) 2*pi*z.^3.*sigma(th,z).*sin(th),0,pi);
else
    disp('dimension d should be either 2 or 3')
end

Sigma = Sigma_hat(zeta)*v/lc; % total scattering cross-section (in 1/sec)

%Sigma = integral(sigma,0,2*pi);
if material.acoustics && d==2
    sigmaNorm = @(th) pi*zeta^2*sigma(th,zeta)/Sigma_hat(zeta);
    pdf = sigmaNorm(xth).*(1/pi);
elseif material.acoustics && d==3
    sigmaNorm = @(th) 4*pi*zeta^3*sigma(th,zeta)/Sigma_hat(zeta);
    pdf = sigmaNorm(xth).*(sin(xth)/2);
else
    disp('dimension d should be either 2 or 3')
end
% sigmaNorm = @(th) sigma(th,zeta)/Sigma;

cdf = cumsum(pdf)*mean(diff(xth));
[cdf,inds] = unique(cdf);
invcdf = griddedInterpolant(cdf,xth(inds));
end


function h = heaviside(x)
h = (x>0);
end