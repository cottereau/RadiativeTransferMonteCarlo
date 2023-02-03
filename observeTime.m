function obs = observeTime(obs,it,P)

% check normalization of dir
r = sqrt(sum(P.x.^2,2));
cospsi = dot(P.x,P.dir,2)./r;
cospsi(cospsi>1)=1;
cospsi(cospsi<-1)=-1;
psi = acos(cospsi);
obs.energy(:,:,it) = obs.energy(:,:,it) + ...
                           histcounts2( r, psi, obs.binX, obs.binPsi )/P.N;
