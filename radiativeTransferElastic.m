function obs = radiativeTransferElastic( source, material, observation )

% constants
T = observation.time(end); % final time of the simulation
Nt = 1000;                 % maximum number of time steps

% initialization
% obs.path: complete paths of some particles
% obs.points: time and direction of passages close to observation sensors
% x,y: current position of each particle
% d: angle of propagation of each particle
% t: current time of each particle
% p: polarization of each particle
obs = initializeObservation(observation,Nt);
material = prepareSigma(material);
[x,y,d,t,p,v,meanFreePath] = initializeParticle(source,material);
it = 1;

% loop on particles
while min(t)<T
    [x,y,t,x0,y0,t0] = propagateParticle(x,y,d,t,v,meanFreePath);
    obs = observePaths(obs,x,y,it,p);
    obs.sensors{it} = observePoints(obs.rp,x,y,x0,y0,t,t0,d,v,p);
    it = it+1;
    if it==Nt, T = min(t); obs.t = obs.t(obs.t<=T); break, end
    [d,p,v,meanFreePath] = scatterParticle(d,material,p,v,meanFreePath);
end

% post-treatment on points
obs = treatmentPoints(obs,source.numberParticles,it);

% checks for debug
if obs.check
    plotPath(obs,source)
end