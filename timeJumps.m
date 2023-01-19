function P = timeJumps(P,dt)

% draw number of jumps for each particle
Nj = poissrnd(dt./P.meanFreePath);
maxNj = max(Nj);

% draw times of jumps
tj = inf(P.N,maxNj);
for i2 = 1:maxNj
    ind = Nj>=i2;
    tj(ind) = rand(sum(ind),1)*dt;
end
tj = sort(tj,2);

% store for particles
P.Nj = Nj;
P.tj = tj;

