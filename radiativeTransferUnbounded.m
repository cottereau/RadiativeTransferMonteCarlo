function obs = radiativeTransferUnbounded( d, source, material, observation )

% physics
acoustics = material.acoustics;

% discretization in packets of particles  (for optimal vectorization)
Npk = 5e4;                             % size of packets (5e4 seems optimal on my computer)
Np = ceil(source.numberParticles/Npk); % number of packets

% OBSERVATION STRUCTURE
% d         : dimension of the problem
% acoustics : true=acoustics, false=elastics
% N         : total number of particles
% t         : time instants
% Nt        : number of time instants
% psi       : propagation directions
% binPsi    : bins for histograms in direction
% Npsi      : number of directions
% binR      : bins for histograms in positions
% r         : sensor positions
% Nr        : number of sensor positions
% dr        : weight of small interval of radius
% energy    : matrix of observations, size [Nr Npsi Nt]
% dE        : energy of a single particle (depends on r)
[ obs, energy, Ec, Npsi, binPsi, Nr, binR, Nt, t ] = ...
                initializeObservation( d, acoustics, observation, Np*Npk );

% prepare scattering cross sections 
material = prepareSigma( material, d );      

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
    % coherent     : false when particle has been scattered at least once
    P = initializeParticle( Npk, d, acoustics, source, material );
    energyi = zeros( Nr, Npsi, Nt, 2 );
    Eci = zeros( Nt, 1 );
    [energyi(:,:,1,:),Eci(1)] = observeTime( d, P, binPsi, binR );

    % loop on time
    for it = 2:Nt

        % propagate particles
        P = propagateParticle( material, P, t(it) );

        % observe energies (as a function of [Psi x t])
        [energyi(:,:,it,:),Eci(it)] = observeTime( d, P, binPsi, binR );

    % end of loop on time
    end
    
    energy = energy + energyi;
    Ec = Ec + Eci;

% end of loop on packages
end
obs.energy = (1./(obs.dr'*obs.N)).*energy;
obs.Ec = Ec;

% energy density as a function of [x t] and [t]
obs.energyDensity = squeeze(sum(obs.energy,[2 4]));
obs.energyDomain = obs.dr*obs.energyDensity;
