function P = propagateParticleSmallDt(mat,boundaries,P,T)

% choice of time sub-step
T = T-mean(P.t);
dt = min(0.01*mat.meanFreeTime(1),T); % with 0.01, we miss scattering events with a rate of 1e-5
Nt = ceil(T/dt);
dt = T/Nt;

% loop on time sub-iterations
for i1 = 1:Nt

    % polarized particles
    p = P.p;
    s = ~P.p;

    % propagate particles
    x1 = P.x;
    P.x(p,:) = P.x(p,:) + (mat.vp.*dt).*P.dir(p,:);
    P.x(s,:) = P.x(s,:) + (mat.vs.*dt).*P.dir(s,:);

    % bounce on boundaries
    for i2 = 1:length(boundaries)
        bnd = boundaries(i2);
        ind = (P.x(:,bnd.dir)-bnd.val).*(bnd.val-x1(:,bnd.dir))>0;
        P.x(ind,bnd.dir) = (2*bnd.val)-P.x(ind,bnd.dir);
        P.dir(ind,bnd.dir) = -P.dir(ind,bnd.dir);
    end

    % choose particles that are scattered at the end of time step
    Nscat = poissrnd(dt/mat.meanFreeTime(1),P.N,1);
    scat = Nscat>0;
    Nscat = sum(scat);

    % draw scattering angle around direction of propagation
    if P.d==2
        phi = randi([0 1],Nscat,1,'logical');
        P.perp(scat,:) = ((-1).^phi).*P.perp(scat,:);
    elseif P.d==3
        phi = (2*pi)*rand(Nscat,1);
        P.perp(scat,:) = cos(phi).*P.perp(scat,:) ...
            + sin(phi).*cross(P.dir(scat,:),P.perp(scat,:));
    end

    % draw scattering angle away from direction of propagation
    if mat.acoustics
        theta = zeros(P.N,1);
        theta(scat) = mat.invcdf(rand(Nscat,1));

    else
        % change polarization
        p = scat & p;
        P.p( p & (rand(P.N,1)>mat.P2P) ) = false;
        s = scat & s;
        P.p( s & (rand(P.N,1)>mat.S2S) ) = true;

        % draw scattering angle away from direction of propagation
        theta = zeros(P.N,1);
        theta( scat ) = rand(Nscat,1);
        theta(  P.p & p ) = mat.invcdf{1,1}( theta(  P.p & p ) );
        theta( ~P.p & p ) = mat.invcdf{1,2}( theta( ~P.p & p ) );
        theta(  P.p & s ) = mat.invcdf{2,1}( theta(  P.p & s ) );
        theta( ~P.p & s ) = mat.invcdf{2,2}( theta( ~P.p & s ) );
    end

    % scattering away from direction of propagation
    dir = P.dir;
    P.dir =  cos(theta).*dir + sin(theta).*P.perp;
    P.dir = P.dir./sqrt(sum(P.dir.^2,2));
    P.perp = -sin(theta).*dir + cos(theta).*P.perp;
    P.perp = P.perp./sqrt(sum(P.perp.^2,2));

    % end of loop on time sub-iterations
end

% % elastics
% else
% 
%     % loop on time sub-iterations
%     for i1 = 1:Nt
% 
%         % polarized particles
%         p = P.p;
%         s = ~P.p;
% 
%         % propagate particles
%         P.x(p,:) = P.x(p,:) + (mat.vp.*dt).*P.dir(p,:);
%         P.x(s,:) = P.x(s,:) + (mat.vs.*dt).*P.dir(s,:);
% 
%         % choose particles that are scattered at the end of time step
%         Nscat = poissrnd(dt/mat.meanFreeTime(1),P.N,1);
%         scat = Nscat>0;
%         Nscat = sum(scat);
% 
%         % draw scattering angle around direction of propagation
%         phi = zeros(P.N,1);
%         if P.d==2
%             phi(scat) = randi([0 1],Nscat,1,'logical');
%             P.perp = ((-1).^phi).*P.perp;
%         elseif P.d==3
%             phi(scat) = (2*pi)*rand(Nscat,1);
%             P.perp = cos(phi).*P.perp + sin(phi).*cross(P.dir,P.perp);
%         end
% 
%         % change polarization
%         p = scat & p;
%         P.p( p & (rand(P.N,1)>mat.P2P) ) = false;
%         s = scat & s;
%         P.p( s & (rand(P.N,1)>mat.S2S) ) = true;
% 
%         % draw scattering angle away from direction of propagation
%         theta = zeros(P.N,1);
%         theta( scat ) = rand(Nscat,1);
%         theta(  P.p & p ) = mat.invcdf{1,1}( theta(  P.p & p ) );
%         theta( ~P.p & p ) = mat.invcdf{1,2}( theta( ~P.p & p ) );
%         theta(  P.p & s ) = mat.invcdf{2,1}( theta(  P.p & s ) );
%         theta( ~P.p & s ) = mat.invcdf{2,2}( theta( ~P.p & s ) );
% 
%         % scattering away from direction of propagation
%         dir = P.dir;
%         P.dir =  cos(theta).*dir + sin(theta).*P.perp;
%         P.dir = P.dir./sqrt(sum(P.dir.^2,2));
%         P.perp = -sin(theta).*dir + cos(theta).*P.perp;
%         P.perp = P.perp./sqrt(sum(P.perp.^2,2));
% 
%     % end of loop on time sub-iterations
%     end
% 
% end

% update particle times
P.t = P.t+T;
