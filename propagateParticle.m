function P = propagateParticle(mat,P,T)

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
    Lj = P.v(ind).*dt(ind);
    dphij = P.dphi(ind);
    dthj = P.dth(ind);
    P.x(ind) = P.x(ind) + Lj.*cos(dphij).*sin(dthj);
    P.y(ind) = P.y(ind) + Lj.*sin(dphij).*sin(dthj);
    P.z(ind) = P.z(ind) + Lj.*cos(dthj);

    % scatter particles (except in last jump)
    P.dphi(ind2) = P.dphi(ind2) + mat.invcdf(rand(sum(ind2),1));
    % this one is commented to stay in 2D -- uncomment when ready
    % P.dth(ind2) = P.dth(ind2) + rand(sum(ind2),1);

    % remaining jumping particles
    P.t(ind) = P.t(ind) + dt(ind);
    dt(ind) = T - P.t(ind);
    ind = dt>0;

end

% compute position in cylindrical coordinates
[P.phi,P.theta,P.r] = cart2sph(P.x,P.y,P.z);
%P.theta = mod(P.theta,pi);
%P.phi = mod(P.phi,2*pi);


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
