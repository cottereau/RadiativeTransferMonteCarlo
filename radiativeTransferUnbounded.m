function obs = radiativeTransferUnbounded( d, source, material, observation )

% physics
acoustics = material.acoustics;

% discretization in packets of particles  (for optimal vectorization)
Npk = 2e4;                             % size of packets (5e4 seems optimal on my computer)
Npk
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
[ obs, Ei, Ec, Npsi, binPsi, Nr, binR, Nt, t ] = ...
                initializeObservation( d, acoustics, observation, Np*Npk );
obs.lambda = source.lambda;
obs.v = material.v;

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
    x = zeros( Npk, 3, Nt );
    dir = zeros( Npk, 3, Nt );
    coherent = false( Npk, Nt );

    % loop on time
    for it = 2:Nt

        % propagate particles
        P = propagateParticle( material, P, t(it) );

        % observe energies (as a function of [Psi x t])
        x(:,:,it) = P.x;
        dir(:,:,it) = P.dir;
        coherent(:,it) = P.coherent;

    % end of loop on time
    end
    
    Ei = Ei + observeTimeNew( d, x, dir, coherent, binPsi, binR );
    Ec = Ec + sum(coherent,1)';

% end of loop on packages
end
obs.energyIncoherent = (1./(obs.dr'*obs.N)).*Ei;
obs.energyDomainCoherent = Ec/obs.N;

% energy density as a function of [x t] and [t]
obs.Ei = squeeze(sum(obs.energyIncoherent,2));
obs.Ec = coherentInABox(obs.energyDomainCoherent,obs.r,0,0,[0 0 0],t,d, ...
                                                 source.lambda,material.v);
obs.energyDomainIncoherent = (obs.dr*obs.Ei)';
