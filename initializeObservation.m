function obs = initializeObservation(observation,material,N)

% basic characteristics of the problem
d = material.dimension;     
acoustics = material.acoustics;

% times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% angles between propagation direction and position vector
binPsi = linspace(0,2*pi,observation.Ndir);
psi = (binPsi(1:end-1)+binPsi(2:end))/2;
Npsi = length(psi);

% sensor positions
binX = observation.sensors;
x = (binX(1:end-1)+binX(2:end))/2;
Nx = length(x);

% initialize matrix of observations
energy = zeros(Nx,Npsi,Nt);

% energy in a small volume of the domain
dV = volumeEnergy(d,x);
dE = 1/N;

% initialize structure
obs = struct('t', t, ...               % time instants
             'Nt', Nt, ...             % number of time instants
             'binPsi', binPsi, ...     % bins for histograms in direction
             'psi', psi, ...           % propagation directions
             'Npsi', Npsi, ...         % number of directions
             'binX', binX, ...         % bins for histograms in positions
             'x', x, ...               % sensor positions
             'Nx', Nx, ...             % number of positions
             'energy', energy, ...     % matrix of observations
             'dV', dV, ...             % small volume of domain
             'dE', dE, ...             % energy of a single particle
             'acoustics', acoustics ); % true=acoustics, false=elastics

end

function dV = volumeEnergy(d,r)
dr = mean(diff(r));
dV = (pi*dr)*(2*r).^(d-1);
end