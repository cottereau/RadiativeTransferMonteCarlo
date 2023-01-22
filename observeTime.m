function obs = observeTime(obs,it,P)

dir = mod(P.d-P.theta,2*pi);
obs.energy(:,:,it) = histcounts2( P.r, dir, obs.binX, obs.binTheta )/P.N;
