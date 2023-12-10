function P = propagateParticle(mat,P,T)

% initialization 
dt = T-P.t;
ind = dt>0;

% loop on flying time until end of interval
while any(ind)

    % number of remaining scattering events for each particle
    Nj = zeros(P.N,1);
    Nj(ind&P.p) = poissrnd(dt(ind&P.p)./mat.meanFreeTime(1));
    if length(mat.meanFreeTime)>1
        Nj(ind&~P.p) = poissrnd(dt(ind&~P.p)./mat.meanFreeTime(2));
    end
    ind2 = Nj>0;
    Nind2 = sum(ind2);
    P.coherent(ind2) = false;

    % flying time until next scattering event
    if any(ind2)
        dt(ind2) = timeNextJump(Nj(ind2),dt(ind2));
    end

    % propagate particles
    p = ind &  P.p;
    P.x(p,:) = P.x(p,:) + (mat.vp.*dt(p)).*P.dir(p,:);
    s = ind & ~P.p;
    P.x(s,:) = P.x(s,:) + (mat.vs.*dt(s)).*P.dir(s,:);

    % scattering angle phi (around direction of propagation)
    if P.d==3
        phi = (2*pi)*rand(Nind2,1);
        P.perp(ind2,:) = cos(phi).*P.perp(ind2,:) ...
                       + sin(phi).*cross(P.dir(ind2,:),P.perp(ind2,:));
    end

    % scattering for acoustics
    if mat.acoustics

        % draw scattering angle theta  (away from direction of propagation)
        theta = mat.invcdf(rand(Nind2,1));

    % scattering for elastics
    else

        % change polarization
        p = ind2 &  P.p;
        P.p( p & (rand(N,1)>mat.P2P) ) = false;
        s = ind2 & ~P.p;
        P.p( s & (rand(N,1)>mat.S2S) ) = true;

        % draw scattering angle theta  (away from direction of propagation)
        theta = rand(N,1);
        theta(  P.p & p ) = mat.invcdf{1,1}( theta(  P.p & p ) );
        theta( ~P.p & p ) = mat.invcdf{1,2}( theta( ~P.p & p ) );
        theta(  P.p & s ) = mat.invcdf{2,1}( theta(  P.p & s ) );
        theta( ~P.p & s ) = mat.invcdf{2,2}( theta( ~P.p & s ) );
 
    end

    % scattering by theta (away from direction of propagation)
    if P.d==2
        rotheta = randi([0 1],Nind2,1,'logical');
        theta(rotheta) = -theta(rotheta);
    end
    P.dir(ind2,:) =  cos(theta).*P.dir(ind2,:) + sin(theta).*P.perp(ind2,:);
    P.perp(ind2,:)= -sin(theta).*P.dir(ind2,:) + cos(theta).*P.perp(ind2,:);

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
