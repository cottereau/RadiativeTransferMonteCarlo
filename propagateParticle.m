function P = propagateParticle(mat,P,T)

% initialization 
dt = T-P.t;
ind = dt>0;

% loop on jumps
while any(ind)

    % select number of jumps on remaining interval dt
    Nj = zeros(P.N,1);
    Nj(ind) = poissrnd(dt(ind)./P.meanFreeTime(ind));
    ind2 = Nj>0;
    Nind2 = sum(ind2);

    % flying time until next jump (or end of interval)
    if any(ind2)
        dt(ind2) = timeNextJump(Nj(ind2),dt(ind2));
    end

    % propagate particles
    P.x(ind,:) = P.x(ind,:) + (P.v(ind).*dt(ind)).*P.dir(ind,:);

    % scatter particles
    theta = mat.invcdf(rand(Nind2,1));
    if P.d==2
        rotheta = randi([0 1],Nind2,1,'logical');
        theta(rotheta) = -theta(rotheta);
    end
    costheta = cos(theta);
    sintheta = sin(theta);
    dir = P.dir(ind2,:);
    perp = P.perp(ind2,:);
    P.dir(ind2,:) = costheta.*dir + sintheta.*perp;
    P.perp(ind2,:) = -sintheta.*dir + costheta.*perp;
    if P.d==3
        phi = 2*pi*rand(Nind2,1);
        P.dir(ind2,:) = rodrigues(P.dir(ind2,:),dir,phi);
        P.perp(ind2,:) = rodrigues(P.perp(ind2,:),dir,phi);
    end
    P.coherent(ind2) = false;

    % remaining jumping particles
    P.t(ind) = P.t(ind) + dt(ind);
    dt(ind) = T - P.t(ind);
    ind = dt>0;

    % plot for debug
    % figure; quiver(P.x(~ind2,1),P.x(~ind2,2),P.dir(~ind2,1),P.dir(~ind2,2),'b')
    % hold on; quiver(P.x(ind2,1),P.x(ind2,2),P.dir(ind2,1),P.dir(ind2,2),'r')

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

% rotation of vector x around axis omega by angle phi (Rodrigues' rotation
% formula)
function y = rodrigues(x,omega,phi)
cosphi = cos(phi);
omegax = dot(omega,x,2).*(1-cos(phi));
y = x.*cosphi + cross(omega,x,2).*sin(phi) + omegax.*omega;
end