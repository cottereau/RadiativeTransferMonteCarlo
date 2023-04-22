function [ obs, energy, Npsi, binPsi, Nr, binR, Nt, t, dE ] = ...
                      initializeObservation( d, acoustics, observation, N )

% times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% angles between propagation direction and position vector
binPsi = linspace(0,pi,observation.Ndir);
psi = (binPsi(1:end-1)+binPsi(2:end))/2;
Npsi = length(psi);

% sensor positions
r = observation.r;
Nr = length(r);
dr = mean(diff(r));
binR = (r(1:end-1)+r(2:end))/2;
binR = [-dr/2 binR binR(end)+dr/2];

% initialize matrix of observations
energy = zeros(Npsi,Nx,Nt);

% energy in a small volume of the domain
[dr,dpsi]= volumeEnergy(d,r,psi);
dE = (1/N)./(dpsi'*dr);

% initialize structure
obs = struct('t', t, ...               % time instants
             'Nt', Nt, ...             % number of time instants
             'binPsi', binPsi, ...     % bins for histograms in direction
             'psi', psi, ...           % propagation directions
             'Npsi', Npsi, ...         % number of directions
             'dpsi', dpsi, ...         % weight of small interval of angle
             'binR', binR, ...         % bins for histograms in positions
             'r', r, ...               % sensor positions
             'Nr', Nr, ...             % number of positions
             'dr', dr, ...             % weight of small interval of radius
             'dE', dE, ...             % energy of a single particle
             'd', d, ...               % dimension of the problem
             'acoustics', acoustics ); % true=acoustics, false=elastics
%             'energy', energy, ...     % matrix of observations

end

function [dr,dphi] = volumeEnergy(d,r,phi)
d0r = mean(diff(r));
dphi = mean(diff(phi));
if d==2
    dr = r*d0r;
    dr(1) = d0r^2/8;
    dphi = 2*pi*ones(size(phi));
elseif d==3
    dr = r.^2*d0r + d0r^2/12;
    dr(1) = d0r^3/24;
    dphi = 2*pi*sin(phi).*dphi;
end
end
