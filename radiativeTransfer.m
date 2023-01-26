function obs = radiativeTransfer( source, material, observation )

% OBSERVATION STRUCTURE
% d         : dimension of the problem
% acoustics : true=acoustics, false=elastics
% t         : time instants
% Nt        : number of time instants
% theta     : propagation directions
% binTheta  : bins for histograms in direction
% Nth       : number of directions
% binX      : bins for histograms in positions
% x         : sensor positions
% Nx        : number of positions
% energy    : matrix of observations size [Nx Nth Nt]
obs = initializeObservation(observation); 
material = prepareSigma(material);        % prepare scattering cross sections 
obs.d = material.dimension;               % space dimensions
obs.acoustics = material.acoustics;       % true for acoustics, false for elastics

% loop on packages of particles
Npk = 5e4;                                % size of packages (5e4 seems optimal on my computer)
Np = ceil(source.numberParticles/Npk);    % number of packages
for ip = 1:Np

    % PARTICLES
    % N            : number of particles
    % x,y          : cartesian coordinates
    % r, theta        : cylindrical coordinates
    % d            : propagation direction (angle between 0 and 2pi)
    % p            : polarization (used only in elasticity)
    % meanFreePath : mean free path
    % v            : propagation velocity
    % Nj           : number of jumps
    % tj           : time of jumps
    P = initializeParticle(Npk,source,material);
    obs = observeTime(obs,1,P);

    % loop on time
    for i1 = 2:obs.Nt

        % draw number of jumps
        dt = obs.t(i1)-obs.t(i1-1);
        P = timeJumps(P,dt);

        % propagate particles
        P = propagateParticle(material,P);

        % observe energies
        obs = observeTime(obs,i1,P);

    % end of loop on time
    end

% end of loop on packages
end

% energy density as a function of [x t]
obs.totalEnergy = squeeze(sum(obs.energy,2));

