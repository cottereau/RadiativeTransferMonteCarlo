function obs = radiativeTransfer( source, material, observation )

% constants
d = material.dimension;                      % space dimensions
acoustics = material.acoustics;              % true for acoustics, false for elastics
T = observation.time(end);                   % final time of the simulation
Nt = 1000;                                   % maximum number of time steps

% pre-processing
% obs.path: complete paths of some particles
% obs.points: time and direction of passages close to observation sensors
obs = initializeObservation(observation,Nt); % preparing output
material = prepareSigma(material);           % prepare scattering cross sections 

% loop on packages of particles
Npk = 5e4;                                   % size of packages (5e4 seems optimal)
Np = ceil(source.numberParticles/Npk);       % number of packages
for ip = 1:Np

    % initialization
    % P current status for each particle:
    %  .x, .y,       : position in cartesian coordinates   % x=r y=0
    %  .r, .theta,   : position in cylindrical coordinates % theta=0
    %  .t,           : time
    %  .d,           : propagation direction (angle between -pi and pi)
    %  .p,           : polarization (used only in elasticity)
    %  .meanFreePath : mean free path
    %  .v            : propagation velocity
    % P0 previous status for each particle
    P = initializeParticle(Npk,source,material);
    obs.sensors = cell(Nt,1);
    it = 1;

    % loop on particles
    while min(P.t)<T
        [P,P0] = propagateParticle(P);
     %   P = position2angle(d,P,false);
     %   obs = observePaths(obs,ip,it,P);
        obs.sensors{it} = observePoints(obs.rp,P,P0);
        it = it+1;
        if it==Nt, T = min(P.t); obs.t = obs.t(obs.t<=T); break, end
        P = scatterParticle(acoustics,material,P);

    end

    % compute energies
    obs = computeEnergy(acoustics,obs,Npk*Np,ip,it);

% end of loop on packages
end

% energy density as a function of [t rp]
obs.totalEnergy = squeeze(sum(obs.energy,1));

% energy on grid
if obs.movie
    obs = computeGridEnergy(obs,acoustics);
end

% % checks for debug
% if obs.check
%     plotPath(obs,source)
% end