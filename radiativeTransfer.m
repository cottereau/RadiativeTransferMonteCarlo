function obs = radiativeTransfer( source, material, observation )

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
obs = initializeObservation(observation); 
obs.d = material.dimension;               % space dimensions
obs.acoustics = material.acoustics;       % true for acoustics, false for elastics
material = prepareSigma(material);        % prepare scattering cross sections 

% loop on packages of particles (for optimal vectorization)
Npk = 5e4;                                % size of packages (5e4 seems optimal on my computer)
Np = ceil(source.numberParticles/Npk);    % number of packages
for ip = 1:Np

    % PARTICLES
    % N            : number of particles
    % x,y,z        : cartesian coordinates
    % r,theta,phi  : cylindrical coordinates
    % d            : direction of propagation
    % p            : polarization (used only in elasticity)
    % meanFreePath : mean free path
    % v            : propagation velocity
    % t            : current time for the particle
    P = initializeParticle(Npk,source,material);
    obs = observeTime(obs,1,P);

    % loop on time
    for i1 = 2:obs.Nt

        % propagate particles
        P = propagateParticle(material,P,obs.t(i1));

        % observe energies (as a function of [x Psi t])
        obs = observeTime(obs,i1,P);

    % end of loop on time
    end

% end of loop on packages
end

% energy density as a function of [x t]
obs.totalEnergy = squeeze(sum(obs.energy,2));

