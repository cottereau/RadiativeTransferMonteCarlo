function E = observeTime(obs,P)

% compute energy normalization
dE = repmat( obs.dE./obs.dV', [1 obs.Npsi]);

% compute radius and angle between direction of propagation and position
% check normalization of dir
r = sqrt(sum(P.x.^2,2));
cospsi = dot(P.x,P.dir,2)./r;
cospsi(cospsi>1)=1;
cospsi(cospsi<-1)=-1;
psi = acos(cospsi);

% accumulate energies
E = histcounts2( r, psi, obs.binX, obs.binPsi ).*dE;
