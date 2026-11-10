function obsWeighted = reweightRadiationPattern(geometry, acoustics, observation, P_history, radPatternFun, modeFractions)
%REWEIGHTRADIATIONPATTERN  Post-process a base MC run to another launch pattern.
%
%   obsWeighted = reweightRadiationPattern(geometry, acoustics, observation, ...
%                                          P_history, radPatternFun, modeFractions)
%
%   Change of measure on the launch direction P.dir0. Within each launched
%   mode beta the particle weight is
%
%       w_i = f_beta * R(dir0_i,beta) / sum_{j in beta} R(dir0_j,beta),
%
%   so each mode keeps the energy fraction f_beta and the weights sum to one.
%   modeFractions = [fP fS] (e.g. src.fP, src.fS of sourcePointForce); the
%   default keeps the base run's launched mode fractions. Densities use the
%   solver's own cell measure (obs.cellMeasure) and resolved axes.
%
%   LIMITATION: a weight cannot change a particle's S polarization. Valid
%   for scalar-S (legacy) transport only; for the polarized solver launch
%   the physical source directly (type 'point_force' / 'moment_tensor').
%
%   P_history{it} is the particle state at observation.time(it); P.dir0 and
%   the launched mode are read from P_history{1} (mode at launch).
%
%   See also SOURCEPOINTFORCE, SOURCEMOMENTTENSOR, OBSERVETIME.

P0 = P_history{1};
Ntotal = P0.N;
launchedP = P0.p(:);
R = max(0, reshape(radPatternFun(P0.dir0, launchedP),[],1));
modes = {launchedP, ~launchedP};
if nargin < 6 || isempty(modeFractions)
    modeFractions = [mean(launchedP) mean(~launchedP)];
end
if acoustics, modes = modes(1); modeFractions = 1; end
weights = zeros(Ntotal,1);
for k = 1:numel(modes)
    in = modes{k};
    if modeFractions(k) == 0, continue; end
    total = sum(R(in));
    if total == 0
        error('Fonte:ZeroPattern','Radiation pattern has zero weight in a mode with positive fraction.');
    end
    weights(in) = modeFractions(k)*R(in)/total;
end
weights = weights/sum(modeFractions);

[obs, E, bins, ibins, vals, Nt] = initializeObservation(geometry, acoustics, observation, Ntotal);
for it = 1:Nt
    P = P_history{it};
    E(:,:,it,:) = observeTime(geometry,acoustics,P.x,P.p,P.dir,bins,ibins,vals,P.alive,weights);
end
obsWeighted = obs;
obsWeighted.weights = weights;
obsWeighted.energyDensity = E./obs.cellMeasure;
obsWeighted.energy = reshape(sum(sum(E,1),2),Nt,[]);
end
