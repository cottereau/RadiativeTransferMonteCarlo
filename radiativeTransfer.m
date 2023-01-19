function obs = radiativeTransfer( source, material, observation )

% constants
T = observation.time(end);                   % final time of the simulation
Nt = 1000;                                   % maximum number of time steps

% pre-processing
obs = initializeObservation(observation,Nt); % preparing output
material = prepareSigma(material);           % prepare scattering cross sections 
obs.d = material.dimension;                  % space dimensions
obs.acoustics = material.acoustics;          % true for acoustics, false for elastics

% loop on packages of particles
Npk = 5e4;                                   % size of packages (5e4 seems optimal)
Np = ceil(source.numberParticles/Npk);       % number of packages
for ip = 1:Np

    % initialization
    % P current status for each particle:
    %  .r            : distance of particle from origin
    %  .d,           : propagation direction (angle between 0 and 2pi)
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
        P = scatterParticle(obs.acoustics,material,P);

    end

    % compute energies
    obs = computeEnergy(obs.acoustics,obs,Npk*Np,it);

% end of loop on packages
end

% energy density as a function of [t rp]
obs.totalEnergy = squeeze(sum(obs.energy,1));
