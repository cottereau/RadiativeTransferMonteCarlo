function P = propagateParticle(mat,P,T)

% initialization 
dt = T-P.t;
ind = dt>0;

% loop on flying time until end of interval
while any(ind)

    % polarized particles
    p = ind &  P.p;
    s = ind & ~P.p;

    % number of remaining scattering events for each particle
    Nj = zeros(P.N,1);
    Nj(p) = poissrnd(dt(p)./mat.meanFreeTime(1));
    if length(mat.meanFreeTime)>1
        Nj(s) = poissrnd(dt(s)./mat.meanFreeTime(2));
    end
    scat = Nj>0;
    Nscatter = sum(scat);
    P.coherent(scat) = false;

    % flying time until next scattering event
    if any(scat)
        dt(scat) = timeNextJump(Nj(scat),dt(scat));
    end

    % propagate particles
    P.x(p,:) = P.x(p,:) + (mat.vp.*dt(p)).*P.dir(p,:);
    P.x(s,:) = P.x(s,:) + (mat.vs.*dt(s)).*P.dir(s,:);

    % scattering around direction of propagation
    if P.d==2
        phi = randi([0 1],Nscatter,1,'logical');
        P.perp(scat,:) = ((-1).^phi).*P.perp(scat,:);
    elseif P.d==3
        phi = (2*pi)*rand(Nscatter,1);
        P.perp(scat,:) = cos(phi).*P.perp(scat,:) ...
                     + sin(phi).*cross(P.dir(scat,:),P.perp(scat,:));
    end

    % scattering for acoustics
    if mat.acoustics

        % draw scattering angle away from direction of propagation
        theta = mat.invcdf(rand(Nscatter,1));

    % scattering for elastics
    else

        % change polarization
        p = scat & p;
        P.p( p & (rand(P.N,1)>mat.P2P) ) = false;
        s = scat & s;
        P.p( s & (rand(P.N,1)>mat.S2S) ) = true;

        % draw scattering angle away from direction of propagation
        theta = rand(P.N,1);
        theta(  P.p & p ) = mat.invcdf{1,1}( theta(  P.p & p ) );
        theta( ~P.p & p ) = mat.invcdf{1,2}( theta( ~P.p & p ) );
        theta(  P.p & s ) = mat.invcdf{2,1}( theta(  P.p & s ) );
        theta( ~P.p & s ) = mat.invcdf{2,2}( theta( ~P.p & s ) );
        theta = theta(scat);

    end

    % scattering away from direction of propagation
    dir = P.dir(scat,:);
    P.dir(scat,:) =  cos(theta).*dir + sin(theta).*P.perp(scat,:);
    P.dir(scat,:) = P.dir(scat,:)./sqrt(sum(P.dir(scat,:).^2,2));
    P.perp(scat,:)= -sin(theta).*dir + cos(theta).*P.perp(scat,:);
    P.perp(scat,:) = P.perp(scat,:)./sqrt(sum(P.perp(scat,:).^2,2));

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
