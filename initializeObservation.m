function obs = initializeObservation(observation,Nt)

% initialization
obs = struct('sensors', [], ...
             'binTime', observation.time, ...
             'binTheta', linspace(0,pi,observation.Ndir+1) );
% vector of times
obs.t = (obs.binTime(1:end-1)+obs.binTime(2:end))/2;

% vector of angle (for direction)
theta = (obs.binTheta(1:end-1)+obs.binTheta(2:end))/2;
obs.theta = [theta pi+theta theta(1)];

% observation at given sensors
if isfield(observation,'sensors')
    rp = observation.sensors;
    if size(rp,1)<size(rp,2); rp = rp'; end
else
    rp = [];
end
obs.rp = rp;

% initializing passages at sensors
% sensors is a cell. Each element corresponds to a time step, and contains
% [radius of sensor, angle of crossing, time of crossing, direction of
% crossing]
obs.sensors = cell(Nt,1);

