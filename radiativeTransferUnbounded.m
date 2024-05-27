function obs = radiativeTransferUnbounded( d, source, material, observation )

% physics
acoustics = material.acoustics;

% discretization in packets of particles  (for optimal vectorization)
Npk = 1e4; % size of packets 
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
[ obs, Ei, Ec, binPsi, binR, Nt, t ] = ...
      initializeObservation( d, acoustics, material, observation, Np*Npk );

% prepare scattering cross sections 
material = prepareSigma( material, d );      

% loop on packages of particles
parfor ip = 1:Np

    % PARTICLES
    % N            : number of particles
    % d            : dimension of the problem
    % x            : cartesian coordinates
    % dir          : direction of propagation
    % perp         : orthogonal to direction of propagation
    % p            : polarization
    % t            : current time for the particle
    % coherent     : false when particle has been scattered at least once
    P = initializeParticle( Npk, d, acoustics, source );
    x = zeros( Npk, 3, Nt ); x(:,:,1) = P.x;
    p = true( Npk, Nt ); p(:,1) = P.p;
    dir = zeros( Npk, 3, Nt ); dir(:,:,1) = P.dir;
    coherent = false( Npk, Nt ); coherent(:,1) = P.coherent;

    % loop on time
    for it = 2:Nt

        % propagate particles
        P = propagateParticle( material, P, t(it) );

        % observe energies (as a function of [Psi x t])
        x(:,:,it) = P.x;
        p(:,it) = P.p;
        dir(:,:,it) = P.dir;
        coherent(:,it) = P.coherent;

    % end of loop on time
    end
    
    Ei = Ei + observeTime( d, acoustics, x, p, dir, coherent, binPsi, binR );
    Ec = Ec + [sum( coherent & p ,1)' sum( coherent & ~p ,1)'];

% end of loop on packages
end
obs.energyIncoherent = (1./(obs.dr'*obs.N)).*double(Ei);
obs.energyDomainCoherent = Ec(:,1:1+~acoustics)/obs.N;

% energy density as a function of [x t] and [t]
obs.Ei = squeeze(sum(obs.energyIncoherent,2));
obs.Ec = coherentInABox(obs.energyDomainCoherent,obs.r,0,0,[0 0 0],t,d, ...
                                                 source.lambda,material);
obs.energyDomainIncoherent = shiftdim(pagemtimes(obs.dr,obs.Ei));
