function E = Hoshiba_RTE_Unbounded_MonteCarlo(t,source_station_dists,source,material,geometry,scat_order)

% This function calculates the energy density of scalar waves propagating
% on a (2D or 3D) random isotropic scattering medium (fullspace), where the 
% total scattering cross-section is a constant (direction-independent) based 
% on Hoshiba 1991 (Simulation of multiple-scattered coda wave excitation based
% on the energy conservation law). 

% Inputs 
% t                    : time vector
% source_station_dists : vector of source-station distances
% source               : source characteristics
% material             : material characteristics
% geometry             : geometry of the propagation domain
% scat_order           : max number of scattering events (default : 30)

% Output
% E                    : energy density at requested source-station distances

if nargin==5
    scat_order = 20; % default value
end

d = geometry.dimension;
v = material.v;
Np = source.numberParticles;
dt = mean(diff(t));
r = source_station_dists;
[m,n] = size(r); if m>n, r=r'; end

Sigma = prepareSigmaOne(material.sigma{1},d);
meanFreePath = 1/Sigma;

% Source position
xs = source.position; 

% Receiver position
xr = xs + [r' zeros(size(r')) zeros(size(r'))];
nb_stations = length(r);

energyDensity = zeros(length(t),nb_stations,scat_order+1);

for k=1:length(r)
    for i=1:scat_order
        s = exprnd(meanFreePath,[Np i]);
        phi = 2*pi*rand(Np,i);
        if d==2
            theta = (pi/2)*ones(Np,i);
            dz = zeros(Np,i);
        elseif d==3
            theta = acos(-1+2*rand(Np,i));
            dz = s.*cos(theta);
        else
            disp('dimension of the problem should be either 2 or 3')
        end
        
        dx = s.*sin(theta).*cos(phi);
        dy = s.*sin(theta).*sin(phi);

        x = xs(1) + cumsum(dx,2);
        y = xs(2) + cumsum(dy,2);
        z = xs(3) + cumsum(dz,2);

        X = sqrt( (x(:,1)-xs(1)).^2 + (y(:,1)-xs(2)).^2 + (z(:,1)-xs(3)).^2 ); 
        Y = sqrt( diff(x,[],2).^2 + diff(y,[],2).^2 + diff(z,[],2).^2 );
        Z = pdist2([x(:,end) y(:,end) z(:,end)],xr(k,:)); 
    
        times = (X+sum(Y,2)+Z)/v;
        if d==2
            D = squeeze(Z);
            probs = exp(-D/meanFreePath)./(2*pi*D);
        elseif d==3
            D = squeeze(Z);
            probs = exp(-D/meanFreePath)./(4*pi*D.^2);
        else
            disp('dimension of the problem should be either 2 or 3')
        end

        for j=1:Np
            ind = (t>times(j)-dt/2 & t<times(j)+dt/2);
            energyDensity(ind,k,i+1) = energyDensity(ind,k,i+1) + probs(j);
        end
    end

    ind0 = (t==r(k)/v);
    if ~isempty(energyDensity(ind0,k,1))
        if d==2
            energyDensity(ind0,k,1) = energyDensity(ind0,k,1) + exp(-v*t(ind0)/meanFreePath)/(2*pi*r(k));
        elseif d==3
            energyDensity(ind0,k,1) = energyDensity(ind0,k,1) + exp(-v*t(ind0)/meanFreePath)/(4*pi*r(k)^2);
        else
            disp('dimension of the problem should be either 2 or 3')
        end
    end

end

E = sum(energyDensity,3)/v/Np/dt;
