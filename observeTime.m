function [E,Ec] = observeTime(d,P,binPsi,binR)

% compute radius and angle between direction of propagation and position
% check normalization of dir
r = vecnorm(P.x,2,2);
cospsi = dot(P.x,P.dir,2)./r;
cospsi(cospsi>1)=1;
cospsi(cospsi<-1)=-1;

% for d==3, the bins correspond to cos(psi)
if d==2
    psi = acos(cospsi);
elseif d==3
    psi = cospsi;
end

% coherent
ind = P.coherent;

% accumulate energies
E = cat( 4, histcounts2( r(ind), psi(ind), binR, binPsi ), ...
            histcounts2( r(~ind), psi(~ind), binR, binPsi ) );
Ec = sum(ind);
