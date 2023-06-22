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
Nr = length(r);
dr = mean(diff(r));

Sigma = prepareSigmaOne(material.sigma); % homogeneous to 1/[L]
meanFreeTime = 1/v/Sigma;
D = v/2/Sigma;

% a (normalized time)      : Sigma*v*t
% b (normalized distance)  : r*Sigma, i.e. source-station distance * scattering cross-section
a = t/meanFreeTime;
b = r*Sigma;
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
    E = Sigma*dr*(b./a)./real(sqrt(1-(b./a).^2)).*exp(real(sqrt(a.^2-b.^2))-a);
    E(heaviside(b-a)) = 0;
    [~,ind] = min((a-b).^2,[],1);
    for i1 = 1:Nt
        E(ind(i1),i1) = E(ind(i1),i1) + exp(-a(i1));
    end
    E(1,1) = 1;
    if nargout==2
        E_diff = (dr*Sigma)*(b./a).*exp(-b.^2./(2*a));
    end
elseif d==3
    G = @(x) exp(x).*sqrt(1+2.026./x);
    E = (4*pi)*b.^2.*real((1-(b./a).^2).^(1/8)).*exp(-a) ...
        .*G(a.*(1-(b./a).^2).^(3/4)).*heaviside(a-b)./(4*pi*a/3).^(3/2);
    aa = repmat(a,[Nr 1]);
    ind = (a==b);
    if ~isempty(E(ind))
        E(ind) = E(ind) + exp(-aa(ind))/Sigma;
    end
    E(:,1) = zeros(Nr,1);
    E(1,1) = 1/Sigma;

    if nargout==2
        E_diff = b.^2.*exp(-d*b.^2./(4*a))./(4*pi*a/d).^(3/2);
    end
else
    disp('Dimension "d" should be either 2 or 3 !')
end

end

function h = heaviside(x)
h = (x>0);
end