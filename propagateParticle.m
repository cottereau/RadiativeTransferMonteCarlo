function [P,P0] = propagateParticle(P)

% initialization
P0 = P;

% length of propagation path (exponentially distributed with coefficient
% meanFreePath
L = -log(1-rand(P.N,1)).*P.meanFreePath;
P.L = L;

% final positions and times of particles
% using Al-Kashi theorem
%P.x = P.x + P.L.*cos(P.d);
%P.y = P.y + P.L.*sin(P.d);
P.t = P.t + L./P.v;
P.r = sqrt(L.^2+P0.r.^2+2*L.*P0.r.*P0.costheta);
P.costheta = (L.^2+P.r.^2-P0.r.^2)./(2*L.*P.r);
P.costheta(P.costheta>1)=1;
P.costheta(P.costheta<-1)=-1;
P.theta = acos(P.costheta);