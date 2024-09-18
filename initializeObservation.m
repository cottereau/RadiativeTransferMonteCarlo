function [ obs, energy, coherent, bins, ibins, vals, Nt, t, d1, d2 ] = ...
            initializeObservation( d, acoustics, observation, N )

% times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% sensor positions : radius
binR = observation.r;
r = (binR(1:end-1)+binR(2:end))/2;
Nr = length(r);
dr = diff(binR);

% sensor positions : angles
% in 2D, theta in [-pi,pi], phi is unused
% in 3D, theta (azimuth) in [-pi,pi], phi (elevation) in [-pi/2,pi/2]
if isfield(observation,'azimuth') && ~isempty(observation.azimuth)
    binTheta = observation.azimuth;
else
    binTheta = [-2*pi 2*pi];
    disp('no azimuth bins defined: integrating over all angles')
end
theta = (binTheta(1:end-1)+binTheta(2:end))/2;
Ntheta = length(theta);
dtheta = diff(binTheta);
binPhi = [-pi/2 pi/2];
if d==3 && isfield(observation,'elevation') && ~isempty(observation.elevation)
    binPhi = observation.elevation;
elseif d==3
    disp('no elevation bins defined: integrating over all angles')
end
phi = (binPhi(1:end-1)+binPhi(2:end))/2;
Nphi = length(phi);
dphi = diff(binPhi);

% angles between propagation direction and position vector
% for d==3, the bins correspond to cos(psi)
if isfield(observation,'directions') && ~isempty(observation.directions)
    binPsi = observation.directions;
else
    binPsi = [0 pi];
    disp('no directions bins defined: integrating over all angles')
end
if  d==3
    binPsi = cos(binPsi(end:-1:1));
end
psi = (binPsi(1:end-1)+binPsi(2:end))/2;
Npsi = length(psi);
dpsi = diff(binPsi);

% energy in a small volume of the domain
[dr,dphi] = volumeEnergy(d,r,dr,phi,dphi);

% bins for histograms have the following dimensions
ibins = find([Nr Ntheta Nphi Npsi]>1);
if isscalar(ibins) && (ibins==1 || ibins==2)
    ibins = [1 2];
end
if length(ibins)>2
    error('histograms can only be constructed along two directions')
end
if all(ibins==[1 2])
    bins = {binR binTheta};
    vals = {binPhi binPsi};
    d1 = dr;
    d2 = dtheta;
elseif all(ibins==[1 3])
    bins = {binR binPhi};
    vals = {binTheta binPsi};
    d1 = dr;
    d2 = dphi;
elseif all(ibins==[1 4])
    bins = {binR binPsi};
    vals = {binTheta binPhi};
    d1 = dr;
    d2 = dpsi;
else
    error('other combinations of bins to be implemented if need be')
end

% initialize matrix of observations
energy = zeros(length(bins{1})-1,length(bins{2})-1,Nt,1+~acoustics,'uint32');
coherent = zeros(1,Nt);

% initialize observation structured array
obs = struct('d', d, ...                 % dimension of the problem
             'acoustics', acoustics, ... % true if acoustics problem
             'N', N, ...                 % total number of particles
             't', t, ...                 % time instants
             'Nt', Nt, ...               % number of time instants
             'r', r, ...                 % sensor positions
             'Nr', Nr, ...               % number of positions
             'dr', dr, ...               % weight of radius bins
             'binR', binR, ...           % radius bins
             'theta', theta, ...         % sensor angle (elevation)
             'Ntheta', Ntheta, ...       % number of elevations
             'dtheta', dtheta, ...       % weight of azimuth bins
             'binTheta', binTheta, ...   % azimuth bins
             'phi', phi, ...             % sensor angle (azimuth)
             'Nphi', Nphi, ...           % number of azimuths
             'dphi', dphi, ...           % weight of elevation bins
             'binPhi', binPhi, ...       % elevation bins
             'psi', psi, ...             % propagation directions
             'Npsi', Npsi, ...           % number of directions
             'dpsi', dpsi, ...           % weight of direction bins
             'binPsi', binPsi );         % direction bins
%             'energy', energy, ...       % matrix of observation
%             'bins', bins, ...           % bins for histograms

end

% volume energy (optimized for efficiency in 3D)
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
