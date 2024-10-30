function [E,E_diff] = analyticalPaasschens(material,observation,geometry)
% Function calculating the solution of the RTE of scalar waves propagating 
% in an isotropic scattering random medium excited by an isotropic (explosion) 
% source with the Diffusion approx. in 2D & 3D based on :
% J. C. J. Paasschens. Solution of the time-dependent Boltzmann equation,
% Phys. Rev. E 56(1), pp. 1135-1141 (1997).

% Note: the solution is exact for 2D and approximate for 3D.

% inputs are similar to those of the main RTE solver --> refer to main.m for details 

% outputs
% Energy densities + Diffusion approximation
% E : energy density integrated over all angles 
    % in 2D : E : 2*pi*r*E*dr
    % in 3D : E : 4*pi*r^2*E*dr
% E_diff : diffusion approximation (same normalization used here)

d = geometry.dimension;
v = material.v;
t = observation.time;
r = observation.x;

if isempty(material.sigma)
    error(['The Differential Scattering Cross-Sections '...
        'has been not defined, please defined it usign DSCS Class'])
end

Sigma = prepareSigmaOne(material.sigma{1},d); % homogeneous to 1/[T]
meanFreeTime = 1/Sigma;

% a (normalized time)      : Sigma*t
a = t/meanFreeTime;
% b (normalized distance)  : r*Sigma/v
b = r*Sigma/v;

[m,n] = size(a);
if m>n, a=a'; end
[m,n] = size(b);
if m<n, b=b'; end

% Gaussian pulse used to define a Dirac delta function
width = 0.0005; % Width of the Gaussian pulse
deltaFunction = @(x,w) exp(-x.^2/(2*w^2))/(sqrt(2*pi)*w);
H = @(x) x>=0; % Modified Heviside function

if d==2
    E = (Sigma/v)^2 .* ...
        (1/2/pi * 1./b.^2 .* exp(-a) .* deltaFunction(a./b-1,width) + ...
         1/2/pi * 1./a .* real((1-(b./a).^2).^(-0.5)) .* exp(-a) .* exp(real(sqrt(a.^2-b.^2))) .* H(a./b-1) );
elseif d==3
    G = @(x) exp(x).*sqrt(1+2.026./x);
    E = (Sigma/v)^3 .* ...
        (1/4/pi * 1./b.^3 .* exp(-a) .* deltaFunction(a./b-1,width) + ...
         real((1-(b./a).^2).^(1/8)).*exp(-a).*G(a.*(1-(b./a).^2).^(3/4))./(4*pi/3*a).^(3/2).*H(a./b-1) );
else
    disp('Dimension "d" should be either 2 or 3 !')
end

if nargout==2
    E_diff = (Sigma/v)^d ./ (4*pi/d*a).^(d/2) .* exp(-d/4*b.^2./a);
end

end