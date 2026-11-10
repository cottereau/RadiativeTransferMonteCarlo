function R = rayleighSpherical(geometry,source,material,observation,opts)
%RAYLEIGHSPHERICAL  Rayleigh-wave channel for an outermost spherical shell.
%
%   R = rayleighSpherical(geometry,source,material,observation,opts)
%   propagates the Rayleigh waves launched by a source sitting on the outer
%   free surface of a spherical shell (or solid sphere) of radius R0, and
%   returns their energy on the grid used by radiativeTransfer with
%   geometry.frame = 'spherical'.
%
%   KINEMATICS: On a sphere of radius R0, Rayleigh waves propagate along
%   great-circle geodesics with constant phase velocity cR. The exact,
%   norm-preserving geodesic update in Cartesian 3D for a time step dt is:
%
%       dtheta = cR * dt / R0;
%       x_new = R0 * ( (x/R0)*cos(dtheta) + v*sin(dtheta) );
%       v_new = -(x/R0)*sin(dtheta) + v*cos(dtheta);
%
%   where x is position (norm R0) and v is unit tangent vector (x . v = 0).
%
%   ENERGY PARTITION AND DEPTH PROFILE:
%   The source partition defaults to Miller & Pursey (1955) on a free
%   surface (67.4% Rayleigh, 25.8% S, 6.9% P for a Poisson solid).
%   The radial depth dependence is computed via RayleighModel, partitioning
%   energy into radial bins according to G(z) = int_0^z (r1^2 + r2^2) dz'.
%
%   opts (all optional):
%     numberParticles  default 2e5
%     partition        [fP fS fR] source energy split, default [.069 .258 .674]
%     substeps         sub-steps per observation time step, default 4
%     direction        'uniform' (isotropic in tangent plane) or 'outgoing'
%
%   See also RAYLEIGHHOMOGENEOUS, RAYLEIGHMODEL, RADIATIVETRANSFER.
%
%   References:
%     - Miller, G.F. & Pursey, H. (1955), On the partition of energy between elastic waves
%       in a semi-infinite solid, Proc. R. Soc. Lond. A 233, 55-69.
%       doi:10.1098/rspa.1955.0245

if nargin < 5, opts = struct; end
def = struct('numberParticles',2e5,'partition',[0.069 0.258 0.674], ...
             'substeps',4,'direction','uniform');
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts,f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
end

%% 1. Geometry: find outermost spherical shell boundary (dir == 5)
R0 = [];
if isfield(geometry,'bnd')
    for i = 1:numel(geometry.bnd)
        b = geometry.bnd(i);
        if b.dir == 5
            if isempty(R0) || b.val > R0
                R0 = b.val;
            end
        end
    end
end
if isempty(R0)
    error('rayleighSpherical expects a spherical boundary with dir = 5.');
end

%% 2. Rayleigh velocity and depth profile
M = RayleighModel(material.vp,material.vs,material.Frequency,0,0,1,R0);
cR = M.cR;
lambdaR = cR / material.Frequency;

% Tabulated depth integral G(z) = int_0^z (r1^2 + r2^2) dz'
Gtab = M.G1 + M.G2;
G = @(z) interp1(M.zgrid, Gtab, min(max(z,0), R0), 'linear');
W = @(d0,d1) max(0, (G(d1) - G(d0)) / Gtab(end)); % energy share in depth interval

%% 3. Observation grid
% Match initializeObservation / observeTime conventions for spherical frame
binR   = observation.x(:).';   % radial bin edges
binTh  = [-pi pi];            % azimuth bin edges (default)
binPhi = [-pi/2 pi/2];        % elevation bin edges (default)

if isfield(observation,'y') && ~isempty(observation.y)
    binTh = observation.y(:).';
end
if isfield(observation,'z') && ~isempty(observation.z)
    binPhi = observation.z(:).';
end

tobs = [0 setdiff(observation.time,0)];
Nt   = numel(tobs);

% Cell centers
rc   = (binR(1:end-1)   + binR(2:end))   / 2;  Nr   = numel(rc);
thc  = (binTh(1:end-1)  + binTh(2:end))  / 2;  Nth  = numel(thc);
phic = (binPhi(1:end-1) + binPhi(2:end)) / 2;  Nphi = numel(phic);

% Volume weights for spherical coordinates: dV = (r_out^3 - r_in^3)/3 * dtheta * (sin(phi_out) - sin(phi_in))
dr3 = (binR(2:end).^3 - binR(1:end-1).^3) / 3;
dth = diff(binTh);
dsinphi = sin(binPhi(2:end)) - sin(binPhi(1:end-1));

% Radial depth weighting for each radial bin [binR(k), binR(k+1)]:
% Depth z is measured inward from outer radius R0:
% z_min = R0 - binR(k+1), z_max = R0 - binR(k)
wDepth = zeros(Nr, 1);
for k = 1:Nr
    z_min = max(0, R0 - binR(k+1));
    z_max = max(0, R0 - binR(k));
    wDepth(k) = W(z_min, z_max);
end
if sum(wDepth) > 0
    wDepth = wDepth / sum(wDepth);
end

%% 4. Source particles on the spherical surface
N  = round(opts.numberParticles);
fR = opts.partition(3);

srcPos = [0 0 R0];
if isfield(source,'position') && ~isempty(source.position)
    srcPos = source.position(:).';
end
if norm(srcPos) == 0
    srcDir = [0 0 1];
else
    srcDir = srcPos / norm(srcPos);
end
srcCenter = R0 * srcDir;

% Build orthonormal tangent basis (e1, e2) at srcCenter
[~, refIdx] = min(abs(srcDir));
refVec = zeros(1,3); refVec(refIdx) = 1;
e1 = cross(srcDir, refVec); e1 = e1 / norm(e1);
e2 = cross(srcDir, e1);     e2 = e2 / norm(e2);

% Launch angles in tangent plane
if strcmpi(opts.direction,'outgoing')
    alpha = zeros(N,1); % along e1
else
    alpha = 2*pi*rand(N,1); % isotropic in tangent plane
end

P_x = repmat(srcCenter, N, 1);
P_v = cos(alpha).*e1 + sin(alpha).*e2;

% If finite source radial profile on surface is given
if isfield(source,'radial') && ~isempty(source.radial)
    rg  = source.radial.GridVectors{1}(:);
    pdf = max(0, source.radial(rg)) .* sin(rg / R0);
    cdf = cumtrapz(rg, pdf); cdf = cdf / cdf(end);
    keep = [true; diff(cdf) > 0];
    dist0 = interp1(cdf(keep), rg(keep), rand(N,1), 'linear', 'extrap');
    % Rotate source positions by dist0 / R0 along random azimuth
    th_rand = 2*pi*rand(N,1);
    tan_rand = cos(th_rand).*e1 + sin(th_rand).*e2;
    ang0 = dist0 / R0;
    P_x = R0 * ( repmat(srcDir, N, 1).*cos(ang0) + tan_rand.*sin(ang0) );
    % Tangent direction update
    P_v = cos(alpha).*e1 + sin(alpha).*e2;
    % Project P_v to be orthogonal to P_x
    n_p = P_x ./ R0;
    P_v = P_v - sum(P_v .* n_p, 2) .* n_p;
    P_v = P_v ./ vecnorm(P_v, 2, 2);
end

%% 5. Geodesic propagation and observation recording
nsub = max(1, round(opts.substeps));

% Output energy matrix matching observation bin dimensions
isAzimuthalIntegrated = (Nth == 1);
if isAzimuthalIntegrated
    E = zeros(Nr, Nphi, Nt);
else
    E = zeros(Nr, Nth, Nphi, Nt);
end

for it = 1:Nt
    if it > 1
        % Advance by this interval, not by mean(diff(tobs)): observation.time
        % is not required to be uniform, and a mean step silently misplaces
        % every particle relative to the time it is recorded at.
        dtheta_step = cR * (tobs(it)-tobs(it-1)) / (nsub * R0);
        for is = 1:nsub
            % Exact geodesic step on 2-sphere
            cos_dt = cos(dtheta_step);
            sin_dt = sin(dtheta_step);
            
            x_norm = P_x ./ R0;
            P_x_new = R0 * ( x_norm .* cos_dt + P_v .* sin_dt );
            P_v_new = -x_norm .* sin_dt + P_v .* cos_dt;
            
            % Numerical drift re-projection
            r_cur = vecnorm(P_x_new, 2, 2);
            P_x = R0 * (P_x_new ./ r_cur);
            n_new = P_x ./ R0;
            P_v = P_v_new - sum(P_v_new .* n_new, 2) .* n_new;
            P_v = P_v ./ vecnorm(P_v, 2, 2);
        end
    end
    
    % Compute spherical coordinates of particles
    [th_p, phi_p, ~] = cart2sph(P_x(:,1), P_x(:,2), P_x(:,3));
    
    % Snap to bin edges to prevent 1-ULP out-of-bounds drops
    phi_p = min(max(phi_p, binPhi(1)), binPhi(end));
    th_p  = min(max(th_p,  binTh(1)),  binTh(end));
    
    % Bin particles in (theta, phi)
    iPhi = discretize(phi_p, binPhi);
    valid = ~isnan(iPhi);
    
    if isAzimuthalIntegrated
        if any(valid)
            cPhi = accumarray(iPhi(valid), 1, [Nphi 1]);
            % Outer product: wDepth (Nr x 1) * cPhi.' (1 x Nphi)
            E(:,:,it) = (wDepth * cPhi.') * (fR / N);
        end
    else
        iTh = discretize(th_p, binTh);
        valid = valid & ~isnan(iTh);
        if any(valid)
            subs = [iTh(valid) iPhi(valid)];
            c2D = accumarray(subs, 1, [Nth Nphi]);
            % 3D accumulation: wDepth(k) * c2D(j, m)
            for k = 1:Nr
                E(k,:,:,it) = wDepth(k) * c2D * (fR / N);
            end
        end
    end
end

%% 6. Compute energy density and format output
R.energy = squeeze(sum(sum(E, 1), 2));
if ~isAzimuthalIntegrated
    R.energy = squeeze(sum(sum(sum(E, 1), 2), 3));
end

% Compute volume-normalized energy density
if isAzimuthalIntegrated
    % dV(k, m) = dr3(k) * diff(binTh) * dsinphi(m)
    dV = (dr3(:) * dsinphi(:).') * diff(binTh);
    R.energyDensity = E ./ reshape(dV, [Nr Nphi 1]);
else
    % 3D dV(k, j, m) = dr3(k) * dth(j) * dsinphi(m)
    dV = zeros(Nr, Nth, Nphi);
    for k = 1:Nr
        for j = 1:Nth
            for m = 1:Nphi
                dV(k,j,m) = dr3(k) * dth(j) * dsinphi(m);
            end
        end
    end
    R.energyDensity = E ./ repmat(dV, [1 1 1 Nt]);
end

R.t = tobs;
R.r = rc;
R.theta = thc;
R.phi = phic;
R.cR = cR;
R.lambdaR = lambdaR;
R.R0 = R0;
R.partition = opts.partition;
R.depthWeights = wDepth;
R.M = M;

end
