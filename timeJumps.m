function P = timeJumps(P,dt)

% draw number of jumps for each particle
Nj = poissrnd(dt./P.meanFreePath);
maxNj = max(Nj);

% draw times of jumps
tj = inf(P.N,maxNj);
for i2 = 1:maxNj
    ind = Nj>=i2;
    tj(ind,i2) = rand(sum(ind),1)*dt;
end
tj = [zeros(P.N,1) dt*ones(P.N,1) tj];
tj = sort(tj,2);

% add to particle structure
P.tj = tj;
P.Nj = Nj;

end