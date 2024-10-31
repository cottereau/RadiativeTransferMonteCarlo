function energyDensity = randomWalkYoshimoto(geometry, source, material, observation, movie)
% This function computes energy densities of acoustic/elastic waves assuming
% isotropic scattering based on a random walk introduced in:
% K. Yoshimoto, Monte Carlo simulation of seismogram envelopes in 
% scattering media. Journal of Geophysical Research: Solid Earth (2000)

% input parameters are the same as those for the radiative transfer solver
% put movie as 'true' to get a video showing the propagating particles

% output
% energyDensity : an array of size [x t] for acoustic and [x t 2] for (P/S) 
%                 elastic waves containing eneregy densities

% Note : this code works only for isotropic scattering

% default: no movie
if nargin<5
    movie = 'false';
end

d = geometry.dimension;
Np = source.numberParticles;   
t = observation.time;        
dt = mean(diff(t));
% sensor positions : radius
binR = observation.r;
r = (binR(1:end-1)+binR(2:end))/2;
dr = diff(binR);

Sigma = prepareSigmaOne(material.sigma{1},d);

if material.acoustics
    
    v = material.v;
    mft = 1/Sigma;           % mean free time
    if dt>0.1*mft
        warning(['Ratio between time step dt and mean free time is large: ...' ...
                  'try to decrease dt'])
    else
        disp(['Ratio between dt and the mean free time is ' num2str(dt/mft)]);
    end

    energyDensity = zeros(length(r),length(t)); % energy density matrix

else

    vp = material.vp; vs = material.vs;
    Sigmap = Sigma(1,1) + Sigma(1,2); Sigmas = Sigma(2,1) + Sigma(2,2);
    mft_PS = 1./sum(Sigma,2);  % mean free times of P & S waves
    
    if dt>0.1*min(mft_PS)
        warning(['Ratio between time step dt and mean free time is large: ...' ...
                  'try to decrease dt'])
    else
        disp(['Ratio between dt and the mean free time is ' num2str(dt/min(mft_PS))]);
    end
    
    % P-to-S & S-to-P conversion probabilities
    pPS = 1-material.P2P;
    pSP = 1-material.S2S;

    energyDensity = zeros(length(r),length(t),2); % energy density matrix

    % Set the mode portions at the source : true = P (explosion), false = S
    if source.polarization == 'P'
        isPwave = true([Np 1]);
    elseif source.polarization == 'S'
        isPwave = false([Np 1]);
    else
        error('Source polarization type should be either "P" or "S"')
    end
end

% Initialize particles (directions & positions)
if d==2
    theta = 2*pi*rand([Np 1]);
    phi = (pi/2)*ones([Np 1]);
elseif d==3
    theta = 2*pi*rand([Np 1]);
    phi = acos(-1+2*rand([Np 1]));
end

% Initial particle positions (random)
radius = abs(randn(Np,1)*source.lambda);
dir = [cos(theta).*sin(phi) sin(theta).*sin(phi) cos(phi)];
xs = radius.*dir;
positions = xs;

if d==2 
    dV = 2*pi*r.*dr;
elseif d==3
    dV = 4*pi*r.^2.*dr;
end

if strcmpi(movie,'true')
    vidfname = 'wave_propagation.avi'; 
    videoObj = VideoWriter(vidfname); % Create video writer object
    open(videoObj);
    fig = figure; % Figure for plotting
end

% Propagate particles
for timeIdx = 1:length(t)
    if material.acoustics
        % Update positions based on wave velocities
        positions = positions + v*dt.*dir;

        % Check for scattering
        scatter_inds = rand(Np,1) < Sigma*dt; % Identify scattering particles
        theta(scatter_inds) = 2*pi*rand([nnz(scatter_inds) 1]);
        if d==3
            phi(scatter_inds) = acos(-1+2*rand([nnz(scatter_inds) 1]));
        end

    else
        % Update positions based on wave type and velocities
        positions = positions + (isPwave .* vp + ~isPwave .* vs) .* dir * dt;
     
        for i = 1:Np
            if isPwave(i) % indicent P wave
                if rand < Sigmap*dt
                    theta(i) = 2*pi*rand; % new direction
                    if d==3
                        phi(i) = acos(-1+2*rand);
                    end
                    if rand < pPS
                        isPwave(i) = false; % convert P to S
                    end
                end
            else % incident S wave
                if rand < Sigmas*dt
                    theta(i) = 2*pi*rand; % new direction
                    if d==3
                        phi(i) = acos(-1+2*rand);
                    end
                    if rand < pSP
                        isPwave(i) = true; % convert S to P
                    end
                end
            end
        end
    end

    % New (or old) propagation direction
    dir = [cos(theta).*sin(phi) sin(theta).*sin(phi) cos(phi)];
    
    % Check if particles are within the receiver volume
    %distances = sqrt(sum((positions-xs).^2,2));
    distances = sqrt(sum(positions.^2, 2));

    if material.acoustics
        % Calculate energy density
        for k=1:length(r)
            withinVolume = r(k)-dr(k)/2 <= distances & distances <= r(k)+dr(k)/2;
            n = sum(withinVolume);
            % Normalization
            energyDensity(k,timeIdx) = n/Np/dV(k);
        end

    else
        % Calculate energy density
        for k=1:length(r)
            withinVolume = r(k)-dr(k)/2 <= distances & distances <= r(k)+dr(k)/2;
            nP = sum(withinVolume & isPwave);
            nS = sum(withinVolume & ~isPwave);
            % Normalization
            energyDensity(k,timeIdx,1) = nP/Np/dV(k); % P wave energy density
            energyDensity(k,timeIdx,2) = nS/Np/dV(k); % S wave energy density
        end
    end

    if strcmpi(movie,'true')
        % Plot particle positions
        clf; hold on; box on;
        if material.acoustics
            scatter(positions(:,1), positions(:,2), 1, 'b', 'filled');
            pbaspect([1 1 1]);
        else
            scatter(positions(isPwave,1), positions(isPwave,2), 1, 'b', 'filled');
            scatter(positions(~isPwave,1), positions(~isPwave,2), 1, 'r', 'filled');
            pbaspect([1 1 1]);
        end
        hold off;
        xlim([-max(r) max(r)]); ylim([-max(r) max(r)]);
        title(['Time: ' num2str(t(timeIdx)) ' s']);
        xlabel('X Position'); ylabel('Y Position');
        if material.acoustics
            legend('Particle');
        else
            legend('P particle','S particle');
        end
        drawnow;
        
        frame = getframe(fig);
        writeVideo(videoObj, frame);
    end
    
    if rem(timeIdx,200)==0
        disp([num2str(timeIdx) ' out of ' num2str(length(t)) ' done!'])
    end
end

if strcmpi(movie,'true')
    close(videoObj);    
    close(fig);
end
