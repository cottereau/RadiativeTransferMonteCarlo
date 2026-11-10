function [ obs, energy, bins, ibins, vals, Nt, t, d1, d2 ] = ...
            initializeObservation( geometry, acoustics, observation, N )

% constants
d = geometry.dimension;
frame = geometry.frame;

% Preserve the requested physical times, including pre-event snapshots.
validateattributes(observation.time,{'numeric'},{'real','finite','vector','nonempty'});
t = double(observation.time(:)');
if any(diff(t)<=0)
    error('RTE:ObservationTimes','Observation times must increase strictly.');
end
Nt = length(t);

% sensor positions : x in cartesian or radius in spherical
binX = observation.x;
binX = validateEdges(binX,'x');
if ~strcmp(frame,'cartesian') && binX(1)<0
    error('RTE:ObservationBins','Radial edges must be nonnegative.');
end
x = (binX(1:end-1)+binX(2:end))/2;
Nx = length(x);
dx = diff(binX);

% sensor positions : (y,z) in cartesian or (theta,phi) in spherical
% in 2D cylindrical, theta in [-pi,pi], phi is unused
% in 3D cylindrical, theta (azimuth) in [-pi,pi], z in [-Inf Inf]
% in 2D spherical, theta in [-pi,pi], phi is unused
% in 3D spherical, theta (azimuth) in [-pi,pi], phi (elevation) in [-pi/2,pi/2]
if strcmp(frame,'cartesian')
    binY = [-Inf Inf];
    binZ = [-Inf Inf];
elseif strcmp(frame,'cylindrical')
    binY = [-pi pi];
    binZ = [-Inf Inf];
elseif strcmp(frame,'spherical')
    binY = [-pi pi];
    binZ = [-pi/2 pi/2];
else
    error('unknown frame type');
end
if isfield(observation,'y') && ~isempty(observation.y)
    binY = observation.y;
else
    disp(['no y bins defined: integrating over [' num2str(binY) ']'])
end
binY = validateEdges(binY,'y');
if ~strcmp(frame,'cartesian') && (any(~isfinite(binY)) || any(abs(binY)>pi))
    error('RTE:ObservationBins','Azimuth edges must lie in [-pi,pi].');
end
y = (binY(1:end-1)+binY(2:end))/2;
Ny = length(y);
dy = diff(binY);
% in 3D spherical, binZ corresponds to sin(phi)
if d==3 && isfield(observation,'z') && ~isempty(observation.z)
    binZ = observation.z;
elseif d==3
    disp(['no z bins defined: integrating over [' num2str(binZ) ']'])
end
if  d==3 && strcmp(frame,'spherical')
    if any(~isfinite(binZ)) || any(abs(binZ)>pi/2)
        error('RTE:ObservationBins','Spherical elevation edges must lie in [-pi/2,pi/2].');
    end
    binZ = sin(binZ);
end
binZ = validateEdges(binZ,'z');
z = (binZ(1:end-1)+binZ(2:end))/2;
Nz = length(z);
dz = diff(binZ);

% angles between propagation direction and position vector
% for d==3, the bins correspond to cos(psi)
if isfield(observation,'directions') && ~isempty(observation.directions)
    binPsi = observation.directions;
else
    binPsi = [0 pi];
    disp(['no directions bins defined: integrating over [' num2str(binPsi) ']'])
end
binPsi = validateEdges(binPsi,'directions');
if any(~isfinite(binPsi)) || binPsi(1)<0 || binPsi(end)>pi
    error('RTE:ObservationBins','Direction edges must lie in [0,pi].');
end
if  d==3
    binPsi = cos(binPsi(end:-1:1));
end
psi = (binPsi(1:end-1)+binPsi(2:end))/2;
Npsi = length(psi);
dpsi = diff(binPsi);

% energy in a small volume of the domain
if strcmp(frame,'spherical') 
    dx = volumeEnergySpherical(d,x,dx);
elseif strcmp(frame,'cylindrical')
    dx = volumeEnergyCylindrical(x,dx);
end

% bins for histograms have the following dimensions
ibins = find([Nx Ny Nz Npsi]>1);

% handle single-bin case (total integration) where ibins is empty
if isempty(ibins)
    % default to tracking X and Y (single bins) to avoid crash
    ibins = [1 2];
elseif isscalar(ibins)
    if ibins~=1
        ibins = [1 ibins];
    else
        ibins = [1 2];
    end
end

if length(ibins) > 2
    error('RTE:ObservationDimensions','The histogram solver supports at most two resolved coordinates.');
end
if strcmp(frame,'spherical') && ~ismember(ibins,[1 2;1 3;1 4],'rows')
    error('RTE:ObservationDimensions','Unsupported spherical coordinate pair.');
end
if strcmp(frame,'cylindrical') && isequal(ibins,[2 3])
    error('RTE:ObservationDimensions','Unsupported cylindrical coordinate pair.');
end

% Use isequal for robust comparison avoiding dimension mismatch errors
if isequal(ibins, [1 2])
    bins = {binX binY};
    vals = {binPsi binZ};
    d1 = dx;
    d2 = dy;
elseif isequal(ibins, [1 3])
    bins = {binX binZ};
    vals = {binY binPsi};
    d1 = dx;
    d2 = dz;
elseif isequal(ibins, [1 4])
    bins = {binX binPsi};
    vals = {binY binZ};
    d1 = dx;
    d2 = dpsi;
elseif isequal(ibins, [2 4])
    bins = {binY binPsi};
    vals = {binX binZ};
    d1 = dy;
    d2 = dpsi;
elseif isequal(ibins, [3 4])
    bins = {binZ binPsi};
    vals = {binX binY};
    d1 = dz;
    d2 = dpsi;
elseif isequal(ibins, [2 3])
    bins = {binY binZ};
    vals = {binX binPsi};
    d1 = dy;
    d2 = dz;

else
    error(['Bin combination [' num2str(ibins) '] not implemented yet']);
end

% initialize matrix of observations
energy = zeros(length(bins{1})-1, length(bins{2})-1, Nt, 1+~acoustics);

% initialize observation structured array
obs = struct('d', d, ...                 % dimension of the problem
             'frame', frame, ...         % observation coordinate frame
             'acoustics', acoustics, ... % true if acoustics problem
             'N', N, ...                 % total number of particles
             't', t, ...                 % time instants
             'Nt', Nt, ...               % number of time instants
             'x', x, ...                 % sensor positions
             'Nx', Nx, ...               % number of positions
             'dx', dx, ...               % weight of radius bins
             'binX', binX, ...           % radius bins
             'y', y, ...                 % sensor angle (elevation)
             'Ny', Ny, ...               % number of elevations
             'dy', dy, ...               % weight of azimuth bins
             'binY', binY, ...           % azimuth bins
             'z', z, ...                 % sensor angle (azimuth)
             'Nz', Nz, ...               % number of azimuths
             'dz', dz, ...               % weight of elevation bins
             'binZ', binZ, ...           % elevation bins
             'psi', psi, ...             % propagation directions
             'Npsi', Npsi, ...           % number of directions
             'dpsi', dpsi, ...           % weight of direction bins
             'binPsi', binPsi );         % direction bins
[obs.cellMeasure,obs.integratedSpatialAxes] = observationCellMeasure(obs,ibins);
obs.resolvedAxes = ibins;
obs.directionResolved = any(ibins==4);
%             'energy', energy, ...       % matrix of observation
%             'bins', bins, ...           % bins for histograms

end

function edges = validateEdges(edges,name)
validateattributes(edges,{'numeric'},{'real','vector','nonempty','nonnan'});
edges = double(edges(:)');
if numel(edges)<2 || any(diff(edges)<=0) || any(isinf(edges(2:end-1)))
    error('RTE:ObservationBins','%s edges must increase strictly with only optional infinite outer edges.',name);
end
end

% volume energy
% total energy in 2D cartesian (the direction variable is psi) 
%       int_x( int_y( int_psi ( a(x,y,psi) dpsi dx dy )))
% total energy in 2D cylindrical/spherical (the direction variable is psi) 
%       int_r( int_th( int_psi ( a(r,th,psi) dpsi dth r dr )))
% total energy in 3D cartesian (the direction variable is cos(psi))
%       int_x( int_y( int_z( int_cospsi ( a(x,y,z,cospsi) dcospsi dth dz dy dx ))))
% total energy in 3D cylindrical (the direction variable is cos(psi))
%       int_r( int_th( int_z( int_cospsi ( a(r,th,z,cospsi) dcospsi dth dz r dr ))))
% total energy in 3D spherical (the direction variable is cos(psi))
%       int_r( int_th( int_sinphi( int_cospsi ( a(r,th,sinphi,cospsi) dcospsi dth dsinphi r^2 dr ))))
function dr = volumeEnergySpherical(d,r,d0r)
if d==2
    dr = r.*d0r;
    if ~isempty(r) && r(1)==0
        dr(1) = d0r(1)^2/8;
    end
else
    dr = r.^2.*d0r + d0r.^3/12;
    if ~isempty(r) && r(1)==0
        dr(1) = d0r(1)^3/24;
    end
end
% if isscalar(r); dr = r^d/d; end
end
function dr = volumeEnergyCylindrical(r,d0r)
dr = r.*d0r;
if ~isempty(r) && r(1)==0
    dr(1) = d0r(1)^2/8;
end
% if isscalar(r); dr = r^2/2; end
end
