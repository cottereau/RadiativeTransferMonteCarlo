function obs = initializeObservation(observation,Nt)

% initialization
obs = struct('sensors', [], ...
             't', observation.time, ...
             'd', linspace(0,pi,observation.Ndir+1), ...
             'movie', observation.movie, ...
             'check', observation.check );

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

% % following some particles
% if obs.check
%     Np = 5;
%     obs.path = zeros(Nt,3,Np);
% end

