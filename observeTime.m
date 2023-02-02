function obs = observeTime(obs,it,P)

dphi = mod(P.dphi-P.phi,2*pi);
obs.energy(:,:,it) = obs.energy(:,:,it) + ...
                histcounts2( P.r, dphi, obs.binX, obs.binPhi )/P.N;
