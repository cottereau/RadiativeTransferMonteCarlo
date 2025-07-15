function plotEnergies( type, obs, material, lambda, cmax, rmax )
% % unbounded case
% if ~isfield(obs,'nSources')
%     obs.nSources = 1;
%     obs.positionSources = [0 0 0];
%     xx = obs.r(obs.r<=max(obs.r)/sqrt(2));
%     obs.boxX = [-xx(end:-1:2) xx];
%     obs.boxZ = obs.boxX;
%     Nx = length(obs.boxX);
%     [boxx,boxz] = meshgrid(obs.boxX,obs.boxZ);
%     r = sqrt(boxx.^2+boxz.^2);
%     E = interp1( obs.r', obs.Ei, r(:), 'linear', 0 );
%     E = E + interp1( obs.r', obs.Ec, r(:), 'linear', 0 );
%     E = permute( reshape(E,Nx,Nx,obs.Nt,Nac), [2 1 3 4] );
%     obs.energyDensityBox = E;
% end

% constants
if nargin<5
    cmax = [];
end

% plot total energy
if isfield(type,'movieTotalEnergy') && type.movieTotalEnergy
    if obs.Npsi==1
        plotTotalEnergy( obs, cmax )
    else
        plotDirectionalEnergy( obs, material, lambda, type.sensors, rmax );
    end
end

% plot total energy and equipartition
if isfield(type,'checkEnergy') && type.checkEnergy
    Etot = sum(obs.energy,3);
    figure; plot(obs.t,Etot,'r-x');
    if ~obs.acoustics
        hold on
        plot(obs.t,obs.energy(:,1)/obs.energy(:,2),'b-')
        eq = (material.vs/material.vp)^(obs.d)/(obs.d-1);
        hold on; plot(obs.t([1 end]),eq*[1 1],'k--')
        legend('Total energy','P/S energy ratio','equipartition ratio (theory)')
    else
        legend('Total energy')
    end
    xlabel('time')
end

if isfield(type,'timehistory') && type.timehistory
    plotTimeHistory(obs,type)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTotalEnergy( obs, cmax )

% values to be plotted
if obs.Nx == 1
    x = obs.y;
    z = obs.z;
else
    x = obs.x;
    if obs.Ny == 1
        z = obs.z;
    else
        z = obs.y;
    end
end
val = obs.energyDensity;

% constants
Nt = length(obs.t);
n = 1+~obs.acoustics;

% estimate cmax and make sure low values are not plotted
if isempty(cmax)
    cmax = squeeze(max(max(obs.energyDensity,[],1),[],2));
    cmax = cmax(ceil(obs.Nt/2),:);
end

% plot total energy - loop on time
figure; lappend = false;
for i1=1:Nt
    ax1 = subplot(n,1,1,'replace');
    surf( ax1, x, z, val(:,:,i1,1)' );
    view(2); shading flat; box on;
    cb1 = colorbar; colormap(ax1,'pink'); clim(ax1,[0 cmax(1)])
    set(ax1,'XLim',[min(x) max(x)],'YLim',[min(z) max(z)],'XTick',[]);
    title(cb1,'P energy')
    title(ax1,['Total Energy Density, time t = ' num2str(obs.t(i1)) 's'])
    if ~obs.acoustics
        ax2 = subplot(n,1,2,'replace');
        surf( ax1, x, z, val(:,:,i1,2)' );
        view(2); shading flat; box on;
        cb2 = colorbar; clim(ax2,[0 cmax(2)]);
        cmap = colormap(ax2,'pink'); colormap(ax2,cmap(end:-1:1,:));
        set(ax2,'XLim',[min(x) max(x)],'YLim',[min(z) max(z)], ...
            'Position',[.13 .23 .6964 .3412]);
        title(cb2,'S energy')
        set(cb2,'position',[.8411 .2298 .05 .31]);
    end
    % export graphics
    exportgraphics(gcf,'movieTotalEnergy.gif','Append',lappend);
    lappend = true;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotDirectionalEnergy( obs, material, lambda, sensors, rmax )

% constants
Ns = size(sensors,1);
Nt = length(obs.t);

% compute directional energy
[psi2pi,Ec,Ei] = directionEnergy( obs, material, lambda, sensors);

% estimate rmax and make sure low values are readable
if nargin<6 | isempty(rmax)
    rmax = squeeze(max(max(max(Ei,[],1),[],2),[],4));
end

% plot directional energy
for i2 = 1:Ns
    % loop on time
    for i1=1:Nt
        if i1==1; figure; lappend = false; else; clf; lappend = true; end
        polarplot(psi2pi,Ei(:,i1,i2,1),'b');
        hold on; polarplot(psi2pi,Ec(:,i1,i2,1),'b--');
        if ~obs.acoustics
            polarplot(psi2pi,Ei(:,i1,i2,2),'r');
            polarplot(psi2pi,Ec(:,i1,i2,2),'r--');
            legend('P - incoherent', 'P - coherent','S - incoherent', 'S - coherent')
        else
            legend('incoherent', 'coherent')
        end
        hold off;
        rlim([0 rmax(i2)])
        title(['sensor at [' num2str(sensors(i2,1)) ',' ...
            num2str(sensors(i2,3)) '], time t = ' num2str(obs.t(i1))])
        legend
        % export graphics
        nameFig = ['movieDirectionalEnergy' num2str(i2) '.gif'];
        exportgraphics(gcf,nameFig,'Append',lappend);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psi2pi,Ec,Ei] = directionEnergy(obs,material,lambda,sensors)
% construct directional energy at sensors

% constant
Ns = size(sensors,1);
[psi,ind] = sort(obs.psi);
psi2pi = [psi 2*pi-psi(end:-1:1)];
if obs.acoustics
    material.vp = material.v;
end
psibounds = [0 (psi2pi(1:end-1)+psi2pi(2:end))/2 2*pi];

% initialization
Ec = zeros(2*obs.Npsi,obs.Nt,Ns,1+~obs.acoustics);
Ei = zeros(2*obs.Npsi,obs.Nt,Ns,1+~obs.acoustics);
obsEc1 = obs.energyCoherent(:,1);
obsEi1 = obs.energyDensityIncoherent(:,ind,:,1);
if ~obs.acoustics
    obsEc2 = obs.energyCoherent(:,2);
    obsEi2 = obs.energyDensityIncoherent(:,ind,:,2);
end

% loop on sources
for i2 = 1:obs.nSources

    % loop on sensors
    for i1 = 1:Ns

        % distance and angle from source to sensor
        X = sensors(i1,:) - obs.positionSources(i2,:);
        r = sqrt(sum(X.^2));
        rxz = sqrt(sum(X([1 3]).^2));
        theta = acos(X(1)/rxz);
        if X(3)<0; theta = 2*pi-theta; end
        psiRot = mod( psi2pi - theta, 2*pi );
        indtheta = theta>=psibounds(1:end-1) & theta<=psibounds(2:end);

        % estimate coherent directional energy at sensor and rotate
        Es = exp(-(r/lambda-(material.vp/lambda)*obs.t).^2 ) .* obsEc1';

        Ec(indtheta,:,i1,1) = Ec(indtheta,:,i1,1) + Es;
        if ~obs.acoustics
            Es = exp(-(r/lambda-(material.vs/lambda)*obs.t).^2 ) .* obsEc2';
            Ec(indtheta,:,i1,2) = Ec(indtheta,:,i1,2) + Es;
        end

        % estimate incoherent directional energy at sensor and rotate
        Es = squeeze( interp1( obs.r', obsEi1, r ));
        Es2 = [ Es(end:-1:1,:,:); Es];
        aux = interp1( psi2pi', Es2, psiRot );
        if any(isnan(aux(:)))
            warning('The directional energy plot was extrapolated, please verify ')
            aux = interp1( psi2pi', Es2, psiRot,'linear','extrap');
        end
        Ei(:,:,i1,1) = Ei(:,:,i1,1) + aux;
        if ~obs.acoustics
            Es = squeeze( interp1( obs.r', obsEi2, r ));
            Es2 = [ Es(end:-1:1,:,:); Es];
            Ei(:,:,i1,2) = Ei(:,:,i1,2) + interp1( psi2pi', Es2, psiRot );
        end

        % end of loop on sources
    end

    % end of loop on sensors
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTimeHistory(obs,type)

% source to sensors distance
s = obs.positionSources(1,:);
di = type.sensors - repmat(s,size(type.sensors,1),1);
r = vecnorm(di')';

n = 1+~obs.acoustics;

% plots
inds = zeros(size(r));
for icap = 1 : numel(r)
    try
        inds(icap) = find(obs.r > r(icap),1,'first');
        a = figure;

        if ~obs.acoustics
            subplot(3,n,1)
            plot( obs.t, obs.Ei(inds(icap),:,1) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Incoherent P Energy')
            subplot(3,2,3)
            plot( obs.t, obs.Ec(inds(icap),:,1) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Coherent P Energy')
            subplot(3,n,5)
            plot( obs.t, obs.Ei(inds(icap),:,1) + obs.Ec(inds(icap),:,1) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Total P Energy')

            subplot(3,n,2)
            plot( obs.t, obs.Ei(inds(icap),:,2) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Incoherent S Energy')
            subplot(3,n,4)
            plot( obs.t, obs.Ec(inds(icap),:,2) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Coherent S Energy')
            subplot(3,n,6)
            plot( obs.t, obs.Ei(inds(icap),:,2) + obs.Ec(inds(icap),:,2) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Total S Energy')

        else
            subplot(3,n,1)
            plot( obs.t, obs.Ei(inds(icap),:) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Incoherent Energy')
            subplot(3,n,2)
            plot( obs.t, obs.Ec(inds(icap),:) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Coherent Energy')
            subplot(3,n,3)
            plot( obs.t, obs.Ei(inds(icap),:) + obs.Ec(inds(icap),:) );
            grid on
            box on
            xlabel('Time [s]')
            ylabel('Total Energy')
        end
        aux = sprintf('Sensor: %g %g %g',(type.sensors(icap,1)),(type.sensors(icap,2)),(type.sensors(icap,3)));
        sgtitle(aux)

        saveas(a,['Sensor_',num2str(icap),'_energy.png'])
    catch
        er = sprintf('The sensor %d has not been found the coords are: [%f %f %f]',icap,type.sensors(icap,:));
        warning(er)
    end
end

end