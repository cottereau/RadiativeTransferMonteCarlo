function obs = initializeObservation(observation)

% vector of times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% vector of directions
binPhi = linspace(0,2*pi,observation.Ndir);
phi = (binPhi(1:end-1)+binPhi(2:end))/2;
Nphi = length(phi);

% vector of positions of sensors
binX = observation.sensors;
x = (binX(1:end-1)+binX(2:end))/2;
Nx = length(x);

% initialize matrix of observations
energy = zeros(Nx,Nphi,Nt);

% initialize structure
obs = struct('t', t, ...               % time instants
             'Nt', Nt, ...             % number of time instants
             'binPhi', binPhi, ...     % bins for histograms in direction
             'phi', phi, ...           % propagation directions
             'Nphi', Nphi, ...         % number of directions
             'binX', binX, ...         % bins for histograms in positions
             'x', x, ...               % sensor positions
             'Nx', Nx, ...             % number of positions
             'energy', energy );       % matrix of observations
             % 'd'                     % dimension of the problem
             % 'acoustics'             % true=acoustics, false=elastics