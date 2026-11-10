function P = propagatePolarizedParticle(mat,P,T)
%PROPAGATEPOLARIZEDPARTICLE Pure-state unraveling of the elastic Born map.
% Unbounded 3D only; mat must be prepared by MaterialClass.prepareSigma.
% Flight times and target modes use the exact rates mat.Sigma. The polar
% angle is drawn from the unpolarized law mat.invcdf{a,b}, the azimuth
% uniformly, and for an S incident polarization e the direction is kept
% with probability e'*R*e/tr(R) <= 1, where tr(R)=2*sigma{a,b}(theta).
% Azimuth-integrated R is isotropic in the transverse plane, so the rate
% into each target mode does not depend on e: the rejection reshapes only
% the angular law. One Kraus outcome l is then drawn with probability
% |K_l e|^2/(e'*R*e); the new S polarization is K_l e/|K_l e|.
% Amplitudes are those of elasticPolarizedKernels, vectorized.
if isempty(mat.Sigma) || ~isequal(size(mat.Sigma),[2 2]) || isempty(mat.invcdf)
    error('Fonte:PolarizedMaterial','Prepare the material with MaterialClass.prepareSigma first.');
end
mix = polarizedKraus(mat);
T = T+zeros(size(P.t));                         % scalar or per-particle target time
Sig = double(mat.Sigma); total = sum(Sig,2); speed = [mat.vp; mat.vs];
remaining = max(0,T-P.t); remaining(~P.alive) = 0;
active = find(remaining>0);
while ~isempty(active)
    a = 2-P.p(active);                          % 1 = P, 2 = S
    flight = -log(rand(numel(active),1))./total(a);
    dt = min(remaining(active),flight);
    P.x(active,:) = P.x(active,:)+(speed(a).*dt).*P.dir(active,:);
    hit = flight<remaining(active);
    remaining(active) = remaining(active)-dt;
    remaining(active(~hit)) = 0;
    ids = active(hit);
    if ~isempty(ids), P = collide(mat,mix,Sig,total,P,ids); end
    active = ids(remaining(ids)>0);
end
P.t(P.alive) = max(P.t(P.alive),T(P.alive));
end

function P = collide(mat,mix,Sig,total,P,ids)
a = 2-P.p(ids);
b = 1+(rand(numel(ids),1).*total(a) > Sig(a,1));
for ia = 1:2
    for ib = 1:2
        g = ids(a==ia & b==ib);
        if isempty(g), continue; end
        n = P.dir(g,:);
        e = P.perp(g,:);
        e = e-sum(e.*n,2).*n; e = e./sqrt(sum(e.^2,2));
        m = zeros(numel(g),3); amp = zeros(numel(g),1+2*(ib==2),3);
        pending = true(numel(g),1);
        while any(pending)
            k = find(pending);
            theta = mat.invcdf{ia,ib}(rand(numel(k),1));
            mk = direction(n(k,:),theta,2*pi*rand(numel(k),1));
            [ak,gain] = polarizedKraus(mat,mix,ia,ib,n(k,:),mk,e(k,:));
            if ia==2
                accept = rand(numel(k),1).*(2*mat.sigma{ia,ib}(theta)) < gain;
            else
                accept = true(numel(k),1); % P incident: gain equals sigma(theta)
            end
            m(k(accept),:) = mk(accept,:);
            amp(k(accept),:,:) = ak(accept,:,:);
            pending(k(accept)) = false;
        end
        power = reshape(sum(amp.^2,2),[],3);
        cumulative = cumsum(power,2);
        pick = 1+sum(rand(numel(g),1).*cumulative(:,3) > cumulative,2);
        pick = min(pick,3);
        P.dir(g,:) = m;
        P.p(g) = ib==1;
        if ib==2
            out = zeros(numel(g),3);
            for l = 1:3
                sel = pick==l; out(sel,:) = amp(sel,:,l);
            end
            normOut = sqrt(sum(out.^2,2));
            [u,~] = frame(m);
            bad = ~(normOut>0);
            out(bad,:) = u(bad,:); normOut(bad) = 1;
            P.perp(g,:) = out./normOut;
        else
            [u,~] = frame(m); P.perp(g,:) = u;
        end
    end
end
end

function m = direction(n,theta,phi)
[u,w] = frame(n);
m = cos(theta).*n+sin(theta).*(cos(phi).*u+sin(phi).*w);
m = m./sqrt(sum(m.^2,2));
end

function [u,w] = frame(n)
[~,j] = min(abs(n),[],2);
ref = zeros(size(n)); ref(sub2ind(size(n),(1:size(n,1))',j)) = 1;
u = cross(n,ref,2); u = u./sqrt(sum(u.^2,2)); w = cross(n,u,2);
end
