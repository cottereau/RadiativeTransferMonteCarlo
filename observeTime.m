function E = observeTime(P,dE,binPsi,binR)

% compute radius and angle between direction of propagation and position
% check normalization of dir
r = sqrt(sum(P.x.^2,2));
cospsi = dot(P.x,P.dir,2)./r;
cospsi(cospsi>1)=1;
cospsi(cospsi<-1)=-1;
psi = acos(cospsi);

% coherent
ind = P.coherent;

% accumulate energies
E = cat( 4, histcounts2( psi(ind), r(ind), binPsi, binR ).*dE, ...
            histcounts2( psi(~ind), r(~ind), binPsi, binR ).*dE );
