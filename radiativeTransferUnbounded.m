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
[ obs, energy, Npsi, binPsi, Nr, binR, Nt, t, dE ] = ...
                initializeObservation( d, acoustics, observation, Np*Npk );

% prepare scattering cross sections 
material = prepareSigma(material);      

% loop on packages of particles
parfor ip = 1:Np

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
    energyi = zeros( Npsi, Nr, Nt );
    energyi(:,:,1) = observeTime( P, dE, binPsi, binR );

    % loop on time
    for it = 2:Nt

        % propagate particles
        P = propagateParticle( material, P, t(it) );

        % observe energies (as a function of [Psi x t])
        energyi(:,:,it) = observeTime( P, dE, binPsi, binR );

    % end of loop on time
    end
    
    energy = energy + energyi;

% end of loop on packages
end
obs.energy = energy;

% energy density as a function of [x t] and [t]
obs.energyDensity = squeeze(tensorprod(obs.dpsi',obs.energy,1));
obs.energyDomain = squeeze(tensorprod(obs.dr',obs.energyDensity,1));
