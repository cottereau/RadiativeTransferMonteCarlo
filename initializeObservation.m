function [ obs, energy, bins, ibins, vals, Nt, t, d1, d2 ] = ...
            initializeObservation( geometry, acoustics, observation, N )

% constants
d = geometry.dimension;

% times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% sensor positions : x in cartesian or radius in spherical
binX = observation.x;
x = (binX(1:end-1)+binX(2:end))/2;
Nx = length(x);
dx = diff(binX);

% sensor positions : (y,z) in cartesian or (theta,phi) in spherical
% in 2D spherical, theta in [-pi,pi], phi is unused
% in 3D spherical, theta (azimuth) in [-pi,pi], phi (elevation) in [-pi/2,pi/2]
if isfield(geometry,'frame') && strcmp(geometry.frame,'cartesian')
    binY = [-Inf Inf];
    binZ = [-Inf Inf];
else
    binY = [-pi pi];
    binZ = [-pi/2 pi/2];
end
if isfield(observation,'y') && ~isempty(observation.y)
    binY = observation.y;
else
    disp(['no y bins defined: integrating over [' num2str(binY) ']'])
end
y = (binY(1:end-1)+binY(2:end))/2;
Ny = length(y);
dy = diff(binY);
if d==3 && isfield(observation,'z') && ~isempty(observation.z)
    binZ = observation.z;
elseif d==3
    disp(['no z bins defined: integrating over [' num2str(binZ) ']'])
end
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
if  d==3
    binPsi = cos(binPsi(end:-1:1));
end
psi = (binPsi(1:end-1)+binPsi(2:end))/2;
Npsi = length(psi);
dpsi = diff(binPsi);

% energy in a small volume of the domain
if ~(isfield(geometry,'frame') && strcmp(geometry.frame,'cartesian'))
    [dx,dz] = volumeEnergy(d,x,dx,z,dz);
end

% bins for histograms have the following dimensions
ibins = find([Nx Ny Nz Npsi]>1);
if isscalar(ibins) && (ibins==1 || ibins==2)
    ibins = [1 2];
end
if length(ibins)>2
    error('histograms can only be constructed along two directions')
end
if all(ibins==[1 2])
    bins = {binX binY};
    vals = {binPsi binZ};
    d1 = dx;
    d2 = dy;
elseif all(ibins==[1 3])
    bins = {binX binZ};
    vals = {binY binPsi};
    d1 = dx;
    d2 = dz;
elseif all(ibins==[1 4])
    bins = {binX binPsi};
    vals = {binY binZ};
    d1 = dx;
    d2 = dpsi;
else
    error('other combinations of bins to be implemented if need be')
end

% initialize matrix of observations
energy = zeros(length(bins{1})-1,length(bins{2})-1,Nt,1+~acoustics,'uint32');

% initialize observation structured array
obs = struct('d', d, ...                 % dimension of the problem
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
%             'energy', energy, ...       % matrix of observation
%             'bins', bins, ...           % bins for histograms

end

% volume energy (optimized for efficiency in 3D) for spherical frame
% total energy in 2D (the direction variable is psi) 
%       int_r( int_psi ( a(r,psi) dpsi ) (r dr) )
% total energy in 3D (the direction variable is cos(psi))
%       int_r( int_cospsi ( a(r,cospsi) dcospsi ) (r^2 dr) )
function [dr,dphi] = volumeEnergy(d,r,d0r,phi,dphi)
if d==2
    dr = r.*d0r;
    if r(1)==0
        dr(1) = d0r(1)^2/8;
    end
elseif d==3
    dr = r.^2.*d0r + d0r.^3/12;
    if r(1)==0
        dr(1) = d0r(1)^3/24;
    end
    dphi = sin(phi+dphi/2)-sin(phi-dphi/2);
end
if isscalar(r); dr = r^d/d; end
if d==2; dphi = 1; end
end
