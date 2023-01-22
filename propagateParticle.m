function P = propagateParticle(mat,P)

% maximum number of jumps
Nj = size(P.tj,2);
dt = max(P.tj(~isinf(P.tj)));

% loop on jumps
for i1 = 2:Nj

    % select jumping particles
    ind = ~isinf(P.tj(:,i1));

    % propagate particles
    tj = P.tj(ind,i1)-P.tj(ind,i1-1);
    Lj = P.v(ind).*tj;
    dj = P.d(ind);
    P.x(ind) = P.x(ind)+Lj.*cos(dj);
    P.y(ind) = P.y(ind)+Lj.*sin(dj);

    % scatter particles (except in last jump)
    ind2 = (P.tj(:,i1)~=dt);
    ind = ind & ind2;
    P.d(ind) = P.d(ind) + mat.invcdf(rand(sum(ind),1));

end

% compute position in cylindrical coordinates
[P.theta,P.r] = cart2pol(P.x,P.y);

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
