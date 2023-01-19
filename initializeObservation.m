function obs = initializeObservation(observation)

% vector of times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% vector of directions
binTheta = linspace(0,2*pi,observation.Ndir);
theta = (binTheta(1:end-1)+binTheta(2:end))/2;
Nth = length(theta);

% vector of positions of sensors
binX = observation.sensors;
x = (binX(1:end-1)+binX(2:end))/2;
Nx = length(x);

% initialize matrix of observations
energy = zeros(Nx,Nth,Nt);

% initialize structure
obs = struct('t', t, ...               % time instants
             'Nt', Nt, ...             % number of time instants
             'binTheta', binTheta, ... % bins for histograms in direction
             'theta', theta, ...       % propagation directions
             'Nth', Nth, ...           % number of directions
             'binX', binX, ...         % bins for histograms in positions
             'x', x, ...               % sensor positions
             'Nx', Nx, ...             % number of positions
             'energy', energy );       % matrix of observations
             % 'd'                     % dimension of the problem
             % 'acoustics'             % true=acoustics, false=elastics