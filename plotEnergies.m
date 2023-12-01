function plotEnergies( obs, v, lambda, cmax, sensors )

% only plot totalEnergy
if nargin<3
    sensors = [];
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
Y = max(obs.energyDensityBox(:))+1;

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
