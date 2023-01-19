function obs = computeEnergy(acoustics,obs,N,it)

% remove any non-necessary time steps
obs.sensors = obs.sensors(1:it-1);

% reshape sensors
obs.sensors = cat(1,obs.sensors{:});

% initialize energies
Np = length(obs.rp);
if ~isfield(obs,'energy')
    obs.energy = zeros(length(obs.binTheta)-1,length(obs.t),Np);
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
               + histcounts2(tmp(:,1),tmp(:,2),obs.binTheta,obs.binTime)/N;
    end
else
    error('elasticity not implemented yet')
end

% erase sensors
obs = rmfield(obs,'sensors');

end
