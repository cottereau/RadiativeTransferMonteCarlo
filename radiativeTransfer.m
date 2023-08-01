function obs = radiativeTransfer( source, material, observation, geometry )

% constants
lambda = source.lambda;
posS0 = source.position;
v = material.v;
type = geometry.type;
siz = geometry.size;

% construct set of virtual sources
Lmax = material.v*max(observation.time);
[ ns, posSources, Rmax ] = virtualSources( geometry, source.position, Lmax);
observation.r = linspace(0,Rmax,ceil(Rmax/observation.dr));
 
% compute solution in full space
obs = radiativeTransferUnbounded( geometry.dimension, source, material, observation );
energyCoherent = squeeze(sum(obs.energy(:,:,:,1),2));
energyIncoherent = squeeze(sum(obs.energy(:,:,:,2),2));
obs.nSources = ns;
obs.positionSources = posSources;

% preparing plotting the energy
dx = ceil( siz(1)/observation.dr );
obs.boxX = linspace( 0, siz(1), dx );
dz = ceil( siz(3)/observation.dr );
obs.boxZ = linspace( -siz(3), 0, dz );
[boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);

% local amplifications at boundaries for slabs and boxes
amp = ones(size(boxx));
if strcmp(type,'halfspace') || strcmp(type,'slab') || strcmp(type,'box')
    amp = amp .* (1 + exp(-(boxz/lambda).^2));
end
if  strcmp(type,'slab') || strcmp(type,'box')
    amp = amp .* (1 + exp(-((boxz+siz(3))/lambda).^2));
end
if  strcmp(type,'box')
    amp = amp .* (1 + exp(-(boxx/lambda).^2));
    amp = amp .* (1 + exp(-((boxx-siz(1))/lambda).^2));
    amp = amp .* (1 + exp(-(posS0(2)/lambda).^2));
    amp = amp .* (1 + exp(-((posS0(2)-siz(3))/lambda).^2));
end
amp = amp(:);

% preparing amplifications of coherent energy for slabs and boxes
theta = zeros(ns);
phi = zeros(ns);
xc = zeros(ns,3);
yc = zeros(ns);
zc = zeros(ns);
tc = zeros(ns);
for i1 = 1:ns
    x = posSources(:,1)-posSources(i1,1);
    y = posSources(:,2)-posSources(i1,2);
    z = posSources(:,3)-posSources(i1,3);
    [theta(:,i1),~] = cart2pol(x,z);
    eta2 = x.*cos(theta(:,i1))+z.*sin(theta(:,i1));
    [phi(:,i1),r] = cart2pol(eta2,y);    
    xc(:,i1) = (posSources(:,1)+posSources(i1,1))/2;
    yc(:,i1) = (posSources(:,2)+posSources(i1,2))/2;
    zc(:,i1) = (posSources(:,3)+posSources(i1,3))/2;
    tc(:,i1) = r/(material.v*2);
end

% construct superposition of total energy for different sources
E = zeros(numel(boxx),obs.Nt);
%maxampC = zeros(ns,obs.Nt);
for i1 = 1:ns
    x = boxx-posSources(i1,1);
    y = posSources(1,2)-posSources(i1,2);
    z = boxz-posSources(i1,3);
    r = sqrt(x.^2+y.^2+z.^2);
    Ecoherent = interp1(obs.r',energyCoherent,r(:),'linear',0);
    Eincoherent = interp1(obs.r',energyIncoherent,r(:),'linear',0);
    ampC = ones(size(Ecoherent));
    for i2 = setdiff(1:ns,i1)
        x = boxx-xc(i2,i1);
        y = posSources(1,2)-yc(i2,i1);
        z = boxz-zc(i2,i1);
        eta2 =  x*cos(theta(i2,i1))+z*sin(theta(i2,i1));
        zeta = -x*sin(theta(i2,i1))+z*cos(theta(i2,i1));
        eta =  eta2*cos(phi(i2,i1))+y*sin(phi(i2,i1));
        xsi = -eta2*sin(phi(i2,i1))+y*cos(phi(i2,i1));
        d2 = zeta.^2+xsi.^2;
        ind = obs.t >= tc(i2,i1);
        ampC(:,ind) = ampC(:,ind) .* (1 + exp(-(eta(:).^2 ...
                +abs(v.^2*(obs.t(ind).^2-tc(i2,i1)^2)-d2(:)))/lambda.^2));
    end
%    maxampC(i1,:) = sum(ampC,1)-4800;
    E = E + Ecoherent.*ampC + Eincoherent.*amp;
end
E = permute(reshape(E,length(obs.boxZ),length(obs.boxX),obs.Nt),[2 1 3]);
obs.energyDensityBox = E;

end
