function [ obs, energy, Ec, binPsi, binR, Nt, t ] = ...
            initializeObservation( d, acoustics, material, observation, N )

% times
t = [0 setdiff(observation.time,0)];
Nt = length(t);

% angles between propagation direction and position vector
% for d==3, the bins correspond to cos(psi)
if d==2
    binPsi = linspace(0,pi,observation.Ndir);
    psi = (binPsi(1:end-1)+binPsi(2:end))/2;
elseif d==3
    binPsi = linspace(-1,1,observation.Ndir);
    psi = acos((binPsi(1:end-1)+binPsi(2:end))/2);
end  
Npsi = length(psi);

% sensor positions
if isfield(observation,'r')
    r = observation.r;
    if r(1)~=0 || any(abs(diff(r)/mean(diff(r))-1)>1e-12)
        error('invalid sensor vector')
    end
else
    if material.acoustics
        Rmax = material.v*max(t);
    else
        Rmax = material.vp*max(t);
    end        
    r = linspace(0,Rmax,ceil(Rmax/observation.dr));
end
Nr = length(r);
dr = mean(diff(r));
binR = (r(1:end-1)+r(2:end))/2;
binR = [-dr/2 binR binR(end)+dr/2];

% initialize matrix of observations
energy = zeros(Nr,Npsi,Nt,1+~acoustics,'uint16');
Ec = zeros(Nt,2);

% energy in a small volume of the domain
dr = volumeEnergy(d,r);

% initialize structure
obs = struct('acoustics', acoustics, ... % true if acoustics problem
             'N', N, ...                 % total number of particles
             't', t, ...                 % time instants
             'Nt', Nt, ...               % number of time instants
             'binPsi', binPsi, ...       % bins for histograms in direction
             'psi', psi, ...             % propagation directions
             'Npsi', Npsi, ...           % number of directions
             'binR', binR, ...           % bins for histograms in positions
             'r', r, ...                 % sensor positions
             'Nr', Nr, ...               % number of positions
             'dr', dr, ...               % weight of small interval of radius
             'd', d);                    % dimension of the problem
%             'energy', energy, ...       % matrix of observations

end

% volume energy (optimized for efficiency in 3D)
% total energy in 2D (the direction variable is psi) 
%       int_r( int_psi ( a(r,psi) dpsi ) (r dr) )
% total energy in 3D (the direction variable is cos(psi))
%       int_r( int_cospsi ( a(r,cospsi) dcospsi ) (r^2 dr) )
function dr = volumeEnergy(d,r)
d0r = mean(diff(r));
if d==2
    dr = r*d0r;
    dr(1) = d0r^2/8;
elseif d==3
    dr = r.^2*d0r + d0r^3/12;
    dr(1) = d0r^3/24;
end
end
