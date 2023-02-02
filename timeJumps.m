function Nj = timeJumps(P,dt)

% draw number of jumps for each particle
Nj = poissrnd(dt./P.meanFreePath);

end