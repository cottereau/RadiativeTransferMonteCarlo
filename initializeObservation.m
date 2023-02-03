function obs = initializeObservation(observation,material,N)

% basic characteristics of the problem
d = material.dimension;     
acoustics = material.acoustics;

% times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% angles between propagation direction and position vector
binPsi = linspace(0,pi,observation.Ndir);
psi = (binPsi(1:end-1)+binPsi(2:end))/2;
Npsi = length(psi);

% sensor positions
binX = observation.sensors;
x = (binX(1:end-1)+binX(2:end))/2;
Nx = length(x);

% initialize matrix of observations
energy = zeros(Npsi,Nx,Nt);

% energy in a small volume of the domain
[dx,dpsi]= volumeEnergy(d,x,psi);
dE = 1/N;

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
             'energy', energy, ...     % matrix of observations
             'dE', dE, ...             % energy of a single particle
             'acoustics', acoustics ); % true=acoustics, false=elastics

end

function [dx,dphi] = volumeEnergy(d,r,phi)
dr = mean(diff(r));
dphi = mean(diff(phi));
if d==2
    dx = r*dr;
    dphi = 2*pi*ones(size(phi));
elseif d==3
    dx = r.^2*dr;
    dphi = 2*pi*sin(phi).*dphi;
end
end