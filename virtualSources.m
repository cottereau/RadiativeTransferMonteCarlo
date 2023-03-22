function [ns,posS,signS,Rmax] = virtualSources( geom, posS, L)

% create box
boxCorners = computBox(geom,posS);

% for full-space problem, only one source
ns = 1;
signS = 1;

% construct family of sources
if isfield( geom, 'planeX' )
    geom.nx = size(geom.planeX,1);
    ind = true;
    i1=1;
    while any(ind)
        S1 = [2*geom.planeX(i1,1)-posS(:,1) posS(:,2:3)];
        S1 = setdiff(S1,posS,'rows');
        ind = distance2box('min',S1,boxCorners)<L;
        posS = [posS; S1(ind,:)];
        % should be a - for dirichlet ???
        % signS = [signS; geom.planeZ(1,2)*signS(ind)];
        signS = [signS; signS(ind)];
        ns = length(signS);
        if size(geom.planeX,1)==1
            break
        end
        i1 = 3-i1;
    end
end
if isfield( geom, 'planeY' )
    geom.ny = size(geom.planeY,1);
end
if isfield( geom, 'planeZ' )
    ind = true;
    i1=1;
    while any(ind)
        S1 = [posS(:,1:2) 2*geom.planeZ(i1,1)-posS(:,3)];
        S1 = setdiff(S1,posS,'rows');
        ind = distance2box('min',S1,boxCorners)<L;
        posS = [posS; S1(ind,:)];
        % should be a - for dirichlet ???
        % signS = [signS; geom.planeZ(1,2)*signS(ind)];
        signS = [signS; signS(ind)];
        ns = length(signS);
        if size(geom.planeZ,1)==1
            break
        end
        i1 = 3-i1;
    end
end

% choose discretization in space
Rmax = distance2box( 'max', posS, boxCorners );
end

function R = distance2box(type,S,C)
switch type
    case 'min'
        ns = size(S,1);
        ind = [1 2; 2 3; 3 4; 4 1];
        R = zeros(ns,4);
        for i1=1:4
            P2 = C(ind(i1,2),:);
            P1 = C(ind(i1,1),:);
            P12 = P2-P1;
            norm12 = sum(P12.^2);
            for i2 = 1:ns
                uh = dot(S(i2,:)-P1,P12)/norm12;
                uh(uh<0) = 0;
                uh(uh>1) = 1;
                R(i2,i1) = vecnorm(P1+uh*P12-S(i2,:));
            end
        end
        R = min(R,[],2);
    case 'max'
        R = 0;
        for i1 = 1:size(S,1)
            for i2 = 1:4
                R = max(R,vecnorm(C(i2,:)-S(i1,:)));
            end
        end
end
end

