function obs = radiativeTransferUnbounded( d, source, material, observation )

% physics
acoustics = material.acoustics;

% discretization in packets of particles  (for optimal vectorization)
Npk = 5e4;                             % size of packets (5e4 seems optimal on my computer)
Np = ceil(source.numberParticles/Npk); % number of packets

% OBSERVATION STRUCTURE
% d         : dimension of the problem
% acoustics : true=acoustics, false=elastics
% t         : time instants
% Nt        : number of time instants
% psi       : propagation directions
% binPsi    : bins for histograms in direction
% Npsi      : number of directions
% binX      : bins for histograms in positions
% x         : sensor positions
% Nx        : number of positions
% energy    : matrix of observations size [Nx Nth Nt]
% dV        : small volume of domain
% dE        : energy of a single particle
obs = initializeObservation( d, acoustics, observation, Np*Npk );
material = prepareSigma(material);        % prepare scattering cross sections 

% loop on packages of particles
for ip = 1:Np

    % PARTICLES
    % N            : number of particles
    % d            : dimension of the problem
    % acoustics    : true=acoustics, false=elastics
    % x            : cartesian coordinates
    % dir          : direction of propagation
    % perp         : orthogonal to direction of propagation
    % p            : polarization (used only in elasticity)
    % meanFreeTime : mean free time
    % v            : propagation velocity
    % t            : current time for the particle
    P = initializeParticle( Npk, d, acoustics, source, material );
    obs.energy(:,:,1) = obs.energy(:,:,1) + observeTime(obs,P);

    % loop on time
    for i1 = 2:obs.Nt

        % propagate particles
        P = propagateParticle(material,P,obs.t(i1));

        % observe energies (as a function of [Psi x t])
        obs.energy(:,:,i1) = obs.energy(:,:,i1) + observeTime(obs,P);

    % end of loop on time
    end

% end of loop on packages
end

% energy density as a function of [x t] and [t]
obs.energyDensity = squeeze(tensorprod(obs.dpsi',obs.energy,1));
obs.energyDomain = squeeze(tensorprod(obs.dx',obs.energyDensity,1));
