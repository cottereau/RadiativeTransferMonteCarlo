function obs = observeTime(obs,it,P)

obs.energy(:,:,it) = histcounts2( P.r, P.d-P.theta, obs.binX, obs.binTheta )/P.N;
