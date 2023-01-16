function obs = computeEnergy(acoustics,obs,N,ip,it)

% discretization in time
dt = diff(obs.t);
xt = [obs.t(1)-dt(1)/2 obs.t(1:end-1)+(dt/2) obs.t(end)+dt(end)/2];
Nt = length(obs.t);

% discretization in space
dd = diff(obs.d);
xd = [obs.d(1)-dd(1)/2 obs.d(1:end-1)+(dd/2) obs.d(end)+dd(end)/2];
Nd = length(obs.d);

% remove any non-necessary time steps
obs.sensors = obs.sensors(1:it-1);
% if ip==1 && obs.check
%     obs.path = obs.path(1:it,:,:);
% end

% reshape sensors
obs.sensors = cat(1,obs.sensors{:});

% initialize energies
Np = length(obs.rp);
if ~isfield(obs,'energy')
    obs.energy = zeros(Nd,Nt,Np);
end

% energy density as a function of [direction t rp]
if acoustics
    for i1 = 1:Np
        tmp = cat(1,obs.sensors{:,i1});
        % remove polarizations
        tmp = tmp(:,1:2);
        % return angles
        tmp(tmp(:,1)>1,1) = 1;
        tmp(tmp(:,1)<-1,1) = -1;
        tmp(:,1) = acos(tmp(:,1));
        obs.energy(:,:,i1) = obs.energy(:,:,i1) ...
                                  + histcounts2(tmp(:,1),tmp(:,2),xd,xt)/N;
    end
else
    error('elasticity not implemented yet')
end

% erase sensors
obs = rmfield(obs,'sensors');

end
