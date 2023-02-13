function P = propagateParticle(mat,P,T)

% constants
d = mat.dimension;

% initialization 
dt = T-P.t;
ind = dt>0;

% loop on jumps
while any(ind)

    % select number of jumps on remaining intervals
    Nj = zeros(P.N,1);
    Nj(ind) = poissrnd(dt(ind)./P.meanFreePath(ind));
    ind2 = Nj>0;

    % flying time until next jump (or end of interval)
    if any(ind2)
        dt(ind2) = timeNextJump(Nj(ind2),dt(ind2));
    end

    % propagate particles
    L = P.v(ind).*dt(ind);
    L = repmat(L,[1 3]);
    P.x(ind,:) = P.x(ind,:) + L.*P.dir(ind,:);

    % scatter particles (except in last jump)
    % select new polarization and angle
    theta = mat.invcdf(rand(sum(ind2),1));
    if d==3
        phi = rand(sum(ind2),1);
    elseif d==2
        phi = zeros(sum(ind2),1);
    end

    % compute new propagation direction
    costheta = repmat(cos(theta),[1 3]);
    sintheta = repmat(sin(theta),[1 3]);
    cosphi = repmat(cos(phi),[1 3]);
    sinphi = repmat(sin(phi),[1 3]);
    dir1 = P.dir(ind2,:);
    dir2 = P.perp(ind2,:);
    dir3 = cross(dir1,dir2);
    perp = cosphi.*dir2+sinphi.*dir3;
    P.dir(ind2,:) = costheta.*dir1 + sintheta.*perp;
    P.perp(ind2,:) = -sintheta.*dir1 + costheta.*perp;

    % remaining jumping particles
    P.t(ind) = P.t(ind) + dt(ind);
    dt(ind) = T - P.t(ind);
    ind = dt>0;

end

end

% draw time of next jump
% if need be, there is a known distribution for the minimum of n uniform
% random variables over [0 1]: pdf=n*(1-x)^(n-1)
function tj = timeNextJump(Nj,dt)
N = size(Nj,1);
m = max(Nj);
tj = inf(N,m);
for i2 = 1:m
    ind = Nj>=i2;
    tj(ind,i2) = rand(sum(ind),1).*dt(ind);
end
tj = sort(tj,2);
tj = tj(:,1);
end

function P = scatterParticle(mat,P)
%function [dtot,p,v,meanFreePath] = scatterParticle(d,mat,p,v,meanFreePath)
N = P.N;
rd = rand(N,1);
th = mat.invcdf(rd);
% % elastics
% else
%     p0 = p;
%     % change polarization
%     probabilityOfChange = mat.P2P*p + mat.S2S*(1-p);
%     change = rand(N,1)>probabilityOfChange;
%     p(change)=~p(change);
%     % velocity for each particle
%     v(change&p) = mat.vp;
%     v(change&~p) = mat.vs;
%     % meanFreePath for each particle
%     meanFreePath(change&p) = mat.meanFreePathP;
%     meanFreePath(change&~p) = mat.meanFreePathS;    
%     % change direction depending on polarizations
%     th = zeros(N,1);
%     th(p0&p) = mat.invcdfPP(rand(nnz(p0&p),1));
%     th(p0&~p) = mat.invcdfPS(rand(nnz(p0&~p),1));
%     th(~p0&p) = mat.invcdfSP(rand(nnz(~p0&p),1));
%     th(~p0&~p) = mat.invcdfSS(rand(nnz(~p0&~p),1));
% end
%P.d = P.d+th;
P.theta = mod(P.theta+th,2*pi);
ind = P.theta>pi;
P.theta(ind) = -(P.theta(ind)-2*pi);
P.costheta = cos(P.theta);
%P.costheta = cos(acos(P.costheta)+th);

end
