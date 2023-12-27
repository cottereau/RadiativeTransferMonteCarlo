function [E,E_diff] = Paasschens_RTE_Unbounded(source, material, observation, geometry)
% Function calculating the solution of the RTE of scalar waves propagating 
% in an isotropic scattering random medium excited by an isotropic (explosion) 
% source with the Diffusion approx. in 2D & 3D based on Paaschens's 1997 paper. 
% Note: the solution is exact for 2D and approximate for 3D.

% inputs
% source, material, observation, geometry: struct arrays --> refer to main.m for details 

% outputs
% Energy densities + Diffusion approximation

d = geometry.dimension;
v = material.v;
t = observation.time;
Nt = length(t);

Lmax = material.v*max(observation.time);
[~,~,Rmax] = virtualSources( geometry, source.position, Lmax);
r = linspace(0,Rmax,ceil(Rmax/observation.dr));
dr = mean(diff(r));

if ~isfield(material,'sigma')
    material.sigma = PSDF2sigma(material,d);
end
Sigma = prepareSigmaOne(material.sigma,d); % homogeneous to 1/[T]
meanFreeTime = 1/Sigma;

% a (normalized time)      : Sigma*t
% b (normalized distance)  : r*Sigma/v
a = t/meanFreeTime;
b = r*Sigma/v;
% E : energy density integrated over all angles 
    % in 2D : E : 2*pi*r*E*dr
    % in 3D : E : 4*pi*r^2*E*dr
% E_diff : diffusion approximation (same normalization used here)

% debug (Regis): I have not modified the E_diff after renormalization

[m,n] = size(a);
if m>n, a=a'; end
[m,n] = size(b);
if m<n, b=b'; end

if d==2
    E = (dr*Sigma/v)*(b./a)./real(sqrt(1-(b./a).^2)).*exp(real(sqrt(a.^2-b.^2))-a);
    [~,ind] = min((a-b).^2,[],1);
    for i1 = 1:Nt
        E(ind(i1),i1) = E(ind(i1),i1) + exp(-a(i1));
    end
    E(1,1) = 1;
    if nargout==2
        E_diff = (dr*Sigma/v)*(b./a).*exp(-b.^2./(2*a));
    end
elseif d==3
    G = @(x) exp(x).*sqrt(1+2.026./x);
    E = (3/2*dr*sqrt(3/pi)*Sigma/v)*sqrt(a).*(b./a).^2.*real((1-(b./a).^2).^(1/8)) ...
        .*exp(-a).*G(a.*(1-(b./a).^2).^(3/4));
    E(heaviside(b-a)) = 0;
    [~,ind] = min((a-b).^2,[],1);
    for i1 = 1:Nt
        E(ind(i1),i1) = E(ind(i1),i1) + exp(-a(i1));
    end
    E(1,1) = 1;
    if nargout==2
        E_diff = (3/2*dr*sqrt(3/pi)*Sigma/v)*sqrt(a).*(b./a).^2.*exp((-3/4)*b.^2./a);
    end
else
    disp('Dimension "d" should be either 2 or 3 !')
end

end

function h = heaviside(x)
h = (x>0);
end