function sensors = observePoints(rp,P,P0)

% initialization
Np = length(rp);
sensors = cell(1,Np);

% loop on sensors positions
for i1 = 1:Np

    % compute intersections
    sensors{i1} = intersectionPathDirect(rp(i1),P0,P);

%     % construct smallest distance from (0,0,0) to the path of each particle,
%     % only on the segment between the initial and final positions
%     % ind returns the indices of particles for which the smallest distance is
%     % neither the initial nor the final positions
%     [th,rh,thetah,~,ind] = isLocalMinimumOfRadiusOnPath(d,P0,P);
% 
%     % find intersections from initial position to intermediate (when applicable)
%     tmp = intersectionPath(rp(i1),P0.r(ind),P0.t(ind),P0.theta(ind),P0.d(ind), ...
%         rh,th,thetah,P.p(ind));
%     P0.r(ind) = rh;
%     P0.t(ind) = th;
%     P0.theta(ind) = thetah;
% 
%     % find intersections for remaining path
%     tmp2 = intersectionPath(rp(i1),P.r,P.t,P.theta,P.d,P0.r,P0.t,P0.theta,P.p);
% 
%     % merging intersections of both parts of the paths
%     sensors{i1} = [tmp; tmp2];

% end of loop on positions
end
end

function sensors = intersectionPathDirect(rp,P0,P)
halfb = P0.r.*P0.costheta;
c = P0.r.^2-rp^2;
quarterDelta = halfb.^2-c;
indDelta = quarterDelta>0;
halfsqrtDelta = sqrt(abs(quarterDelta));
L1 = -halfb-halfsqrtDelta;
sensors1 = sensorsSolution1(P,P0.r,rp,L1,indDelta);
L2 = -halfb+halfsqrtDelta;
sensors2 = sensorsSolution1(P,P0.r,rp,L2,indDelta);
sensors = [sensors1;sensors2];
end

function sensors = sensorsSolution1(P,r0,rp,L1,indDelta)
ind = indDelta & (L1>0) & (L1<=P.L);
ll = L1(ind);
costheta = (ll+(rp^2-r0(ind).^2)./ll)./(2*rp);
sensors = [costheta P.t(ind) P.p(ind)];
end

function sensors = intersectionPathDirect2(rp,P0,P)

% % case of P on sensor
% indon = abs(P.r-rp)<1e-8; if any(indon); sum(indon), end
% sensors0 = [P.d(indon)-P.theta(indon) P.t(indon) P.p(indon)];
% 
% using Al-Kashi's theorem
% constants
r2 = P.r.^2;
h2 = (P.x-P0.x).^2+(P.y-P0.y).^2;

% finding roots stating intersections
b = -1-(r2-P0.r.^2)./h2;
c = (r2-rp^2)./h2;
D2 = b.^2-4*c;
indreal = D2>0;
sqrtDd2 = zeros(size(D2));
sqrtDd2(indreal) = sqrt(D2(indreal))/2;
%sqrtDd2 = real(sqrt(D2)/2);

% first root: check whether on path
alpha = -b/2+sqrtDd2;
ind = indreal & alpha>=0 & alpha<=1;
% ind = ~indon & indreal & alpha>=0 & alpha<=1;
sensors1 = checkonPath(rp,h2(ind),alpha(ind),P.t(ind),P.r(ind),P.d(ind), ...
    P.theta(ind),P.p(ind,1),P0.t(ind),P0.r(ind),P0.theta(ind));
% if sum(ind)==Inf
%     xp = rp*cos(sensors1(:,1));
%     yp = rp*sin(sensors1(:,1));
%     x = P.r(ind).*cos(P.theta(ind));
%     y = P.r(ind).*sin(P.theta(ind));
%     x0 = P0.r(ind).*cos(P0.theta(ind));
%     y0 = P0.r(ind).*sin(P0.theta(ind));
%     figure;hold on; for i1=1:sum(ind)
%         plot([x0(i1) x(i1)],[y0(i1) y(i1)]);
%     end
%     scatter([x0;xp;x],[y0;yp;y],50,'r');
% end
% second root: check whether on path
alpha = -b/2-sqrtDd2;
%ind =  ~indon & indreal & alpha>=0 & alpha<=1;
ind = indreal & alpha>=0 & alpha<=1;
sensors2 = checkonPath(rp,h2(ind),alpha(ind),P.t(ind),P.r(ind),P.d(ind), ...
    P.theta(ind),P.p(ind,1),P0.t(ind),P0.r(ind),P0.theta(ind));

% merge sensors
sensors = [sensors1;sensors2];
%sensors = [sensors0;sensors1;sensors2];
end

function sensors = checkonPath(rp,h2,alpha,t,r,d,theta,p,t0,r0,theta0)
tr = t+(t0-t).*alpha;
thetar = acos(max(min((r.^2+rp^2-alpha.^2.*h2)./((2*rp)*r),1),-1));
% if ~isreal(thetar)
%     disp('ok')
% end
ind = theta>theta0;
thetar(ind) = theta(ind)-thetar(ind);
thetar(~ind) = theta(~ind)+thetar(~ind);
%indcos = thetar>max(theta,theta0) | thetar<min(theta,theta0);
%thetar(indcos) = -thetar(indcos);
% sensors = [thetar tr p];
% if size(sensors,1)>0
%     disp('ok')
%     xp = rp*cos(thetar);
%     yp = rp*sin(thetar);
%     x = r.*cos(theta);
%     y = r.*sin(theta);
%     x0 = r0.*cos(theta0);
%     y0 = r0.*sin(theta0);
%     figure;hold on; for i1=1:5
%         plot([x0(i1) x(i1)],[y0(i1) y(i1)]);
%         scatter([x0(i1);xp(i1);x(i1)],[y0(i1);yp(i1);y(i1)],50,'r');
%     end
% end
sensors = [d-thetar tr p];
end

% find intersection along a path, assuming the relation between radius and
% coordinate along the path is a bijection
function sensors = intersectionPath(rp,r0,t0,theta0,d0,r,t,theta,p)
sensors = zeros(0,3);
alpha = (rp-r0)./(r-r0);
ind = (alpha>=0 & alpha<=1);
alphaind = alpha(ind);
t0ind = t0(ind);
t0ind = t0ind +alphaind.*(t(ind)-t0ind);
d0ind = d0(ind);
theta0ind = theta0(ind);
% attention: est-ce que cette relation linéaire entre angle et rayon
% est correcte ?
theta0ind = theta0ind +alphaind.*(theta(ind)-theta0ind);
sensors = [sensors; [d0ind-theta0ind t0ind p(ind,1)]];
end

function [th,rh,thetah,phih,ind] = isLocalMinimumOfRadiusOnPath(d,P0,P)
if d==2
    b = - [P.x P.y];
    a = [P0.x P0.y] + b;
elseif d==3
    b = - [P.x P.y P.z];
    a = [P0.x P0.y P0.z] + b;
end
L2 = sum(a.^2,2);
alpha = sum(a.*b,2)./L2;
ind = alpha>0 & alpha<1;
aind = a(ind,:);
alphaind = alpha(ind);
xh = P.x(ind)+alphaind.*aind(:,1);
yh = P.y(ind)+alphaind.*aind(:,2);
if d==2
     [thetah,rh] = cart2pol(xh,yh);
     phih = [];
elseif d==3
     zh = P.z(ind)+alphaind.*aind(:,3);
     [thetah,phih,rh] = cart2sph(xh,yh,zh);
end
Ptind = P.t(ind);
th = Ptind+alphaind.*(P0.t(ind)-Ptind);
end
