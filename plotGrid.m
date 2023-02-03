function M = plotGrid(type,obs,cmax)
% energy on grid
if ~isfield(obs,'gridEnergy')
    obs = computeGridEnergy(obs,obs.acoustics);
end
figure;
Nt = length(obs.t);
switch type
    case 'half'
        v = VideoWriter('halfSpace.avi');
        open(v);
        ind = obs.z<=0 & obs.z>=obs.zmin;
        x = obs.x;
        z = obs.z(ind);
        M(Nt) = struct('cdata',[],'colormap',[]);
        for i1=1:Nt
            surf(obs.x,z,obs.energyHS(:,ind,i1)');
            view(2); colorbar; shading flat
            set(gca,'xlim',[min(x) max(x)],'ylim',[min(z) max(z)])
            set(gca,'PlotBoxAspectRatio',[range(obs.x) range(obs.z(ind)) 1])
            title(['time T = ' num2str(obs.t(i1)) 's'])
            clim([0 cmax])
            M(i1) = getframe;
            writeVideo(v,M(i1));
        end
        close(v);
    case 'full'
        v = VideoWriter('fullSpace.avi');
        open(v);
        x = obs.grid;
        mM = [min(x) max(x)];
        M(Nt) = struct('cdata',[],'colormap',[]);
        for i1=1:Nt
            surf(x,x,obs.gridEnergy(:,:,i1));
            view(2); colorbar; shading flat
            set(gca,'xlim',mM,'ylim',mM)
            set(gca,'PlotBoxAspectRatio',[range(x) range(x) 1])
            title(['time T = ' num2str(obs.t(i1)) 's'])
            clim([0 cmax])
            M(i1) = getframe;
            writeVideo(v,M(i1));
        end
        close(v);
    otherwise
        error('unknown type of problem')
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