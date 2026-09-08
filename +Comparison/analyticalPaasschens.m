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
t = observation.time(:).';
r = ((observation.x(1:end-1)+observation.x(2:end))/2).';

% dissipation
Q = Inf;
if isprop(material,'Q') && ~isempty(material.Q)
    Qval = material.Q(:);
    Q = Qval(1);
end

omega = 2*pi*material.Frequency;
if isfinite(Q)
    absorptionFactor = exp(-omega*t/Q);   % row vector, 1 x Nt
else
    absorptionFactor = ones(size(t));
end

% Note : Equivalent Paasschens absorption length:
% la = v*Q/omega
% absorptionFactor = exp(-v*t/la)

if isempty(material.sigma)
    error(['The Differential Scattering Cross-Sections '...
           'has been not defined, please defined it usign DSCS class'])
end

Sigma = MaterialClass.prepareSigmaOne(material.sigma{1},d); % homogeneous to 1/[T]
meanFreeTime = 1/Sigma;

% a (normalized time)      : Sigma*t
a = t/meanFreeTime;
% b (normalized distance)  : r*Sigma/v
b = r*Sigma/v;

% At t = 0, the point-source solution is a Dirac distribution and has no
% finite value at the radial bin centers. Evaluate the formulas only for
% strictly positive times and leave the corresponding t = 0 column at zero.
positiveTime = t > 0;
aPositive = a(positiveTime);
E = zeros(length(r),length(t));

% Gaussian pulse used to define a Dirac delta function
width = 0.0005; % Width of the Gaussian pulse
deltaFunction = @(x,w) exp(-x.^2/(2*w^2))/(sqrt(2*pi)*w);
H = @(x) x>=0; % Modified Heviside function

if d==2
    if any(positiveTime)
        E(:,positiveTime) = (Sigma/v)^2 .* ...
            (1/2/pi * 1./b.^2 .* exp(-aPositive) .* deltaFunction(aPositive./b-1,width) + ...
             1/2/pi * 1./aPositive .* real((1-(b./aPositive).^2).^(-0.5)) .* ...
             exp(-aPositive) .* exp(real(sqrt(aPositive.^2-b.^2))) .* ...
             H(aPositive./b-1) );
    end
elseif d==3
    G = @(x) exp(x).*sqrt(1+2.026./x);
    if any(positiveTime)
        E(:,positiveTime) = (Sigma/v)^3 .* ...
            (1/4/pi * 1./b.^3 .* exp(-aPositive) .* deltaFunction(aPositive./b-1,width) + ...
             real((1-(b./aPositive).^2).^(1/8)).*exp(-aPositive).* ...
             G(aPositive.*(1-(b./aPositive).^2).^(3/4))./ ...
             (4*pi/3*aPositive).^(3/2).*H(aPositive./b-1) );
    end
else
    disp('Dimension "d" should be either 2 or 3 !')
end

if nargout > 1
    E_diff = zeros(size(E));
    if any(positiveTime)
        E_diff(:,positiveTime) = (Sigma/v)^d ./ ...
            (4*pi/d*aPositive).^(d/2) .* ...
            exp(-d/4*b.^2./aPositive);
    end
end

E = E .* absorptionFactor;
if nargout > 1
    E_diff = E_diff .* absorptionFactor;
end
