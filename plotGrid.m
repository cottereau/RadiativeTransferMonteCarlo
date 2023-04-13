function plotGrid(obs,cmax)
% energy on grid
if ~isfield(obs,'gridEnergy')
    obs = computeGridEnergy(obs,obs.acoustics);
end
figure;
Nt = length(obs.t);
x = obs.boxX;
z = obs.boxZ;
%axis tight manual
%set(gca,'nextplot','replacechildren');
for i1=1:Nt
    surf(x,z,obs.energyDensityBox(:,:,i1)');
    view(2); colorbar; shading flat
    set(gca,'xlim',[min(x) max(x)],'ylim',[min(z) max(z)])
    set(gca,'PlotBoxAspectRatio',[range(x) range(z) 1])
    title(['time T = ' num2str(obs.t(i1)) 's'])
    clim([0 cmax])
    exportgraphics(gcf,'testAnimated.gif','Append',true);
end
end

function obs = computeGridEnergy(obs,acoustics)
Np = length(obs.x);
Nt = length(obs.t);
L = max(obs.x)/sqrt(2);
xg = linspace(0,L,min(4*Np,100));
xg = [-xg(end:-1:2) xg];
Ng = length(xg);
[X,Y] = ndgrid(xg,xg);
rg = sqrt(X.^2+Y.^2);
obs.grid = xg;
obs.gridEnergy = zeros(Ng,Ng,Nt);
for i1 = 1:Nt
    obs.gridEnergy(:,:,i1) = interp1(obs.x,obs.energyDensity(:,i1),rg);
end
if ~acoustics
    obs.gridEnergyP = zeros(length(xg),length(xg),Nt);
    for i1 = 1:Nt
        obs.gridEnergyP(:,:,i1) = interp1(obs.x,obs.sensorEnergyP(:,i1),rg);
    end
    obs.gridEnergyS = zeros(length(xg),length(xg),Nt);
    for i1 = 1:Nt
        obs.gridEnergyS(:,:,i1) = interp1(obs.x,obs.sensorEnergyS(:,i1),rg);
    end
end
end