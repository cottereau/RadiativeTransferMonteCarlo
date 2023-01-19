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
        for i1=1:Nt
            surf(obs.x,z,obs.energyHS(:,ind,i1)');
            view(2); colorbar; shading flat
            set(gca,'xlim',[min(x) max(x)],'ylim',[min(z) max(z)])
            set(gca,'PlotBoxAspectRatio',[range(obs.x) range(obs.z(ind)) 1])
            title(['time T = ' num2str(obs.t(i1)) 's'])
            caxis([0 cmax])
            M(i1) = getframe;
            writeVideo(v,M(i1));
        end
        close(v);
    case 'full'
        v = VideoWriter('fullSpace.avi');
        open(v);
        x = obs.grid;
        mM = [min(x) max(x)];
        for i1=1:Nt
            surf(x,x,obs.gridEnergy(:,:,i1));
            view(2); colorbar; shading flat
            set(gca,'xlim',mM,'ylim',mM)
            set(gca,'PlotBoxAspectRatio',[range(x) range(x) 1])
            title(['time T = ' num2str(obs.t(i1)) 's'])
            caxis([0 cmax])
            M(i1) = getframe;
            writeVideo(v,M(i1));
        end
        close(v);
    otherwise
        error('unknown type of problem')
end
