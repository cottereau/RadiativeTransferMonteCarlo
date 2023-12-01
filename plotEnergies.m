function plotEnergies( obs, v, lambda, cmax, sensors )

% only plot totalEnergy
if nargin<5
    sensors = [];
end

% unbounded case
if ~isfield(obs,'nSources')
    obs.nSources = 1;
    obs.positionSources = [0 0 0];
    xx = obs.r(obs.r<=max(obs.r)/sqrt(2));
    obs.boxX = [-xx(end:-1:2) xx];
    obs.boxZ = obs.boxX;
    Nx = length(obs.boxX);
    [boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);
    r = sqrt(boxx.^2+boxz.^2);
    Ei = interp1(obs.r',obs.Ei,r(:),'linear',0);
    Ec = coherentInABox(obs.energyDomainCoherent,boxx(:),0,boxz(:), ...
                                             [0 0 0],obs.t,obs.d,lambda,v);
    E = permute(reshape(Ec+Ei,Nx,Nx,obs.Nt),[2 1 3]);
    obs.energyDensityBox = E;
end

% compute directional energy
[psi2pi,Ec,Ei] = directionEnergy( obs, v, lambda, sensors );

% constants
Ns = size(sensors,1);
Nt = length(obs.t);
x = obs.boxX;
z = obs.boxZ;
lappend = false;
rmax = 1e-2;
numTile = [4 7+(1:(Ns-1))];
Y = max(Ec(:)+Ei(:))+1;

% initialize figure
til = tiledlayout(2+ceil(max(0,Ns-2)./4),3+(Ns>0));
titleS = cell(Ns,1);
for i1 = 1:Ns
    titleS{i1} = ['sensor at [' num2str(sensors(i1,1)) ',' ...
                                               num2str(sensors(i1,3)) ']'];
end

% loop on time
for i1=1:Nt

    % plot total energy
    if i1==1; nexttile([2 3]); else; nexttile(1); end
    hold off; surf(x,z,obs.energyDensityBox(:,:,i1)');
    view(2); colorbar; shading flat
    set(gca,'xlim',[min(x) max(x)],'ylim',[min(z) max(z)])
    set(gca,'PlotBoxAspectRatio',[range(x) range(z) 1])
    clim([0 cmax])
    title('total energy density')
    if Ns>0
        hold on; scatter3(sensors(:,1),sensors(:,3),Y,50,'r','filled');
    end
    hold on; scatter3(obs.positionSources(1,1),obs.positionSources(1,3),Y,50,'k','filled');

    % loop on sensors
    for i2 = 1:Ns
        nexttile(numTile(i2))
        polarplot(psi2pi,Ei(:,i1,i2),'b');
        hold on; polarplot(psi2pi,Ec(:,i1,i2),'r'); hold off;
        rlim([0 rmax])
        title(titleS{i2})
    end

    % export graphics
    title(til,['time t = ' num2str(obs.t(i1)) 's'])
    exportgraphics(gcf,'movieEnergy.gif','Append',lappend);
    lappend = true;
end

end
