function energyDensity = randomWalkYoshimoto_beta(geometry, source, material, observation, movie)
% This function computes energy densities of acoustic/elastic waves assuming
% isotropic scatter_indsing based on a random walk introduced in:
% K. Yoshimoto, Monte Carlo simulation of seismogram envelopes in
% scatter_indsing media. Journal of Geophysical Research: Solid Earth (2000)

% input parameters are the same as those for the radiative transfer solver
% put movie as 'true' to get a video showing the propagating particles

% output
% energyDensity : an array of size [x t] for acoustic and [x t 2] for (P/S)
%                 elastic waves containing eneregy densities

% Note : this code works only for isotropic scatter_indsing

% default: no movie
if nargin<5
    movie = 'false';
end

d = geometry.dimension;
Np = source.numberParticles;
t = observation.time;
dt = mean(diff(t));
if isfield(observation,'r')
    binR = observation.r;
else
    binR = observation.x;
end
r = (binR(1:end-1)+binR(2:end))/2;  % sensor positions : radius

if ~isempty(material.Sigma)
    Sigma = material.Sigma;
else
    Sigma = prepareSigmaOne(material.sigma{:},d);
end

invcdf = scatteringParameters(material.sigma{1},d);

if material.acoustics
    v = material.v;
    mft = 1/Sigma;  % mean free time
    if dt>0.1*mft
        warning(['Ratio between time step dt and mean free time is large: ...' ...
            'try to decrease dt'])
    else
        disp(['Ratio between dt and the mean free time is ' num2str(dt/mft)]);
    end

    energyDensity = zeros(length(r),length(t)); % energy density matrix
else
    vp = material.vp; vs = material.vs;
    Sigmap = sum(Sigma(1,:)); Sigmas = sum(Sigma(2,:));
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

% Initialize particles (directions)
if d==2
    theta = 2*pi*rand(Np,1); % Angle in [0, 2π)
    dir = [cos(theta) sin(theta)]; % Direction vectors in 2D
elseif d==3
    phi = 2*pi*rand(Np,1);
    theta = acos(-1+2*rand(Np,1));
    dir = [sin(theta).*cos(phi) sin(theta).*sin(phi) cos(theta)];
end

% Initialize particles (positions)
if ~isfield(source,'radial') || isempty(source.radial)
    radius = abs(randn(Np,1)*source.lambda/2);
else
    Rmax = source.radial.GridVectors{1}(end);
    invcdfsource = inverseCDF(source.radial,d,Rmax);
    radius = invcdfsource(rand(N,1));
end

% Default source position
if ~isfield(source,'position')
        source.position = zeros(1,d);
end

if d==2
    if numel(source.position) ~= 2
        error('For this 2D problem, the source position should have 2 components');
    end
    positions = source.position + radius.*[cos(theta) sin(theta)];
elseif d==3
    if numel(source.position) ~= 3
        error('For this 3D problem, the source position should have 3 components');
    end
    positions = source.position + ...
        radius.*[sin(theta).*cos(phi) sin(theta).*sin(phi) cos(theta)];
end

% initial positions of the particles (source)
initialPositions = positions;

if d==2
    dV = pi*(binR(2:end).^2 - binR(1:end-1).^2);
elseif d==3
    dV = (4*pi/3)*(binR(2:end).^3 - binR(1:end-1).^3);
end

if strcmpi(movie,'true')
    vidfname = 'wave_propagation.avi';
    videoObj = VideoWriter(vidfname); % Create video writer object
    open(videoObj);
    fig = figure; % Figure for plotting
end

% Propagate particles
for timeIdx = 1:length(t)

    % Calculate energy density/densities(P/S)
    if material.acoustics

        positions = positions + v*dt.*dir; % update positions based on wave velocities

        % Check for scatter_indsing
        scatter_inds = rand(Np,1) < Sigma*dt; % Identify scattered particles
        Nscatter = sum(scatter_inds);

        % Update directions for scattered particles
        if Nscatter > 0
            dir(scatter_inds,:) = updateDirections(dir(scatter_inds,:), invcdf, d);
        end

        % Calculate distances from initial positions
        distances = sqrt(sum((positions-initialPositions).^2, 2));
        
        % Calculate energy density
        n = histcounts(distances, binR); % counts of particles in each spatial bin
        energyDensity(:, timeIdx) = n(:)/Np./dV(:);  % update energy density
    else

        positions = positions + (isPwave.*vp + ~isPwave.*vs).*dir*dt; % update positions

        % Scattering probabilities for each particle
        scattering_probability = zeros(Np, 1);
        scattering_probability(isPwave)  = Sigmap*dt;  % P-wave scattering probability
        scattering_probability(~isPwave) = Sigmas*dt; % S-wave scattering probability

        scatter_inds = rand(Np, 1) < scattering_probability;
        Nscatter = sum(scatter_inds);

        % Update directions for scattered particles
        if Nscatter > 0
            % Update directions for scattering particles
            dir(scatter_inds, :) = updateDirections(dir(scatter_inds, :), invcdf, d);

            % Mode conversion for scattering particles
            isPwave_scatter = isPwave(scatter_inds);
            mode_conversion_probability = zeros(Nscatter, 1);
            mode_conversion_probability(isPwave_scatter) = pPS;
            mode_conversion_probability(~isPwave_scatter) = pSP;
            mode_conversion = rand(Nscatter, 1) < mode_conversion_probability;
            isPwave(scatter_inds) = xor(isPwave_scatter, mode_conversion);
        end

        % Calculate distances from initial positions
        distances = sqrt(sum((positions-initialPositions).^2, 2));

        % Compute histogram counts for P- and S-wave particles
        nP = histcounts(distances(isPwave), binR);  % counts of P-wave particles in each bin
        nS = histcounts(distances(~isPwave), binR); % counts of S-wave particles in each bin

        % Normalize and update energy density
        energyDensity(:, timeIdx, 1) = nP(:)/Np./dV(:); % P wave energy density
        energyDensity(:, timeIdx, 2) = nS(:)/Np./dV(:); % S wave energy density
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

end

function invcdf = scatteringParameters(sigma, d)
Nth = 1e6;
if d==2
    xth = linspace(0, 2*pi, Nth);
    Sigma = integral(sigma, 0, 2*pi);
    sigmaNorm = @(th) (1/Sigma)*sigma(th);
elseif d==3
    xth = linspace(0, pi, Nth);
    Sigma = 2*pi*integral(@(th) sigma(th).*sin(th), 0, pi);
    sigmaNorm = @(th) (2*pi/Sigma)*sin(th).*sigma(th);
end
pdf = sigmaNorm(xth);
cdf = cumsum(pdf)*(pi/Nth);
ind = find(diff(cdf)>1e-12);
ind = unique([ind ind+1]);
invcdf = griddedInterpolant(cdf(ind),xth(ind));
end


% compute the cumulative distribution radial function corresponding to a
% given (positive) function to draw randomly from it
function invcdf = inverseCDF(f,d,Rmax)
Nth = 10000;
xth = linspace(0,Rmax,Nth);
intF = trapz(xth,(xth.^(d-1)).*f(xth));
pdf = (xth.^(d-1)).*f(xth)/intF;
cdf = cumsum(pdf)*mean(diff(xth));
cdf(1) = 0;
cdf(end) = 1;
ind = find(diff(cdf)>0);
ind = unique([ind ind+1]);
invcdf = griddedInterpolant(cdf(ind),xth(ind),'linear','nearest');
end

function dir_new = updateDirections(dir_old, invcdf, d)
% Function to update the directions of scattered particles
% dir_old: Nxd array of current directions for N particles in d dimensions
% invcdf: inverse cumulative distribution function for scattering angles
% d: dimension (2 or 3)
% Returns:
% dir_new: Nxd array of updated directions

Nscatter = size(dir_old, 1);

if d==2
    % Angle of rotation
    delta_theta = invcdf(rand(Nscatter,1)); % sample scattering angles
    % Current directions
    dir_x = dir_old(:,1);
    dir_y = dir_old(:,2);
    % Rotate directions
    dir_new(:,1) = cos(delta_theta).*dir_x - sin(delta_theta).*dir_y;
    dir_new(:,2) = sin(delta_theta).*dir_x + cos(delta_theta).*dir_y;
elseif d == 3
    % Generate random scattering angles in the local coordinate frame
    delta_theta = invcdf(rand(Nscatter, 1));
    delta_phi = 2*pi*rand(Nscatter, 1);

    % Get two perpendicular vectors to dir_old to form a local basis
    % You can use the Gram-Schmidt process or simply construct them manually
    perp1 = zeros(Nscatter, 3);
    perp2 = zeros(Nscatter, 3);

    for i = 1:Nscatter
        % Get the current direction
        d_old = dir_old(i,:);

        % Generate an arbitrary vector that is not collinear with d_old
        arbitrary = [1, 0, 0];
        if abs(dot(arbitrary, d_old)) > 0.9 % If too aligned, choose another
            arbitrary = [0, 1, 0];
        end

        % Apply Gram-Schmidt to find perp1
        perp1(i,:) = arbitrary - dot(arbitrary, d_old) * d_old;
        perp1(i,:) = perp1(i,:) / norm(perp1(i,:)); % Normalize

        % perp2 is simply the cross product of d_old and perp1
        perp2(i,:) = cross(d_old, perp1(i,:));
        perp2(i,:) = perp2(i,:)/norm(perp2(i,:));
    end

    % Calculate the new direction in the local coordinate system
    dir_new = zeros(Nscatter, 3);
    for i = 1:Nscatter
        % Spherical coordinates in the local frame
        dir_new(i,:) = sin(delta_theta(i)) * cos(delta_phi(i)) * perp1(i,:) + ...
                       sin(delta_theta(i)) * sin(delta_phi(i)) * perp2(i,:) + ...
                       cos(delta_theta(i)) * dir_old(i,:);
    end
    % Normalize dir_new to ensure unit vectors
    dir_new = dir_new ./ sqrt(sum(dir_new.^2, 2));
else
    error('Dimension d must be 2 or 3');
end

end