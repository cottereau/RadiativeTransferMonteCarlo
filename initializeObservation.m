function [ obs, energy, Npsi, binPsi, Nx, binX, Nt, t, dE ] = ...
                      initializeObservation( d, acoustics, observation, N )

% times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% angles between propagation direction and position vector
binPsi = linspace(0,pi,observation.Ndir);
psi = (binPsi(1:end-1)+binPsi(2:end))/2;
Npsi = length(psi);

% sensor positions
x = observation.sensors;
Nx = length(x);
dx = mean(diff(x));
binX = (x(1:end-1)+x(2:end))/2;
binX = [-dx/2 binX binX(end)+dx/2];

% binX = observation.sensors;
% x = (binX(1:end-1)+binX(2:end))/2;

% initialize matrix of observations
energy = zeros(Npsi,Nx,Nt);

% energy in a small volume of the domain
[dx,dpsi]= volumeEnergy(d,x,psi);
dE = (1/N)./(dpsi'*dx);

% initialize structure
obs = struct('t', t, ...               % time instants
             'Nt', Nt, ...             % number of time instants
             'binPsi', binPsi, ...     % bins for histograms in direction
             'psi', psi, ...           % propagation directions
             'Npsi', Npsi, ...         % number of directions
             'dpsi', dpsi, ...         % weight of small interval of angle
             'binX', binX, ...         % bins for histograms in positions
             'x', x, ...               % sensor positions
             'Nx', Nx, ...             % number of positions
             'dx', dx, ...             % weight of small interval of radius
             'dE', dE, ...             % energy of a single particle
             'd', d, ...               % dimension of the problem
             'acoustics', acoustics ); % true=acoustics, false=elastics
%             'energy', energy, ...     % matrix of observations

end

function [dx,dphi] = volumeEnergy(d,r,phi)
dr = mean(diff(r));
dphi = mean(diff(phi));
if d==2
    dx = r*dr;
    dx(1) = dr^2/8;
    dphi = 2*pi*ones(size(phi));
elseif d==3
    dx = r.^2*dr + dr^2/12;
    dx(1) = dr^3/24;
    dphi = 2*pi*sin(phi).*dphi;
end
end
