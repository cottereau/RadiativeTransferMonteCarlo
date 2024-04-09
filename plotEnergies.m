function plotEnergies( type, obs, material, lambda, cmax, rmax )

% constants
Nac = size(obs.energyDomainCoherent,2);

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
    E = interp1( obs.r', obs.Ei, r(:), 'linear', 0 );
    E = E + interp1( obs.r', obs.Ec, r(:), 'linear', 0 );
    E = permute( reshape(E,Nx,Nx,obs.Nt,Nac), [2 1 3 4] );
    obs.energyDensityBox = E;
end


% plot total energy
if isfield(type,'movieTotalEnergy') && type.movieTotalEnergy
    if nargin<5; cmax = []; end
    plotTotalEnergy( obs, cmax );
end

% plot directional energy
if isfield(type,'movieDirectionalEnergy') && type.movieDirectionalEnergy
    if nargin>3 && ~isempty(type.sensors)
        if nargin<6; rmax = []; end
        plotDirectionalEnergy( obs, material, lambda, type.sensors, rmax );
    end
end

% plot total energy and equipartition
if isfield(type,'equipartition') && type.equipartition
    E = obs.energyDomainCoherent+obs.energyDomainIncoherent;
    figure; plot(obs.t,E/max(sum(E,2)));
    if ~obs.acoustics
        eq = (material.vs/material.vp)^(obs.d)/(obs.d-1);
        hold on; plot(obs.t,sum(E,2)/max(sum(E,2)),'k--')
        hold on; plot(obs.t,eq*(E(:,2)./E(:,1)))
        legend('P energy','S energy','total energy','equipartition')
    end
    if ~obs.acoustics
        % Determine the polarization type of the source
        if obs.energyDomainCoherent(1,1)==1 && obs.energyDomainCoherent(1,2)==0
            pol_type = 'P';
        elseif obs.energyDomainCoherent(1,1)==0 && obs.energyDomainCoherent(1,2)==1
            pol_type = 'S';
        else
            error('The source polarization should be either P or S')
        end
        title(['Energy partitioning for a source with polarization type "' pol_type '"'])
        set(gca,'YLim',[0 1.2]);
    else
        title('Total energy')
    end
    xlabel('time')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTotalEnergy( obs, cmax )

% constants
Nt = length(obs.t);
x = obs.boxX;
z = obs.boxZ;
n = 1+~obs.acoustics;
indp = z>=0;
indn = z<0;

if obs.acoustics
    indp = z==z;
end

% estimate cmax and make sure low values are not plotted
if isempty(cmax)
    cmax = squeeze(max(max(obs.energyDensityBox,[],1),[],2));
    cmax = cmax(ceil(obs.Nt/2),:);
end

% plot total energy - loop on time
figure; lappend = false;
for i1=1:Nt
    ax1 = subplot(n,1,1,'replace');
    surf(ax1, x,z(indp),obs.energyDensityBox(:,indp,i1,1)');
    view(2); shading flat; box on;
    cb1 = colorbar; colormap(ax1,'pink'); clim(ax1,[0 cmax(1)])
    set(ax1,'XLim',[min(x) max(x)],'YLim',[min(z(indp)) max(z(indp))], ...
            'XTick',[],'Position',[.13 .5838 .6964 .3412]);
    title(cb1,'P energy')
    set(cb1,'position',[.8411 .5833 .05 .31]);
    title(ax1,['Total Energy Density, time t = ' num2str(obs.t(i1)) 's'])
    if ~obs.acoustics
        ax2 = subplot(n,1,2,'replace');
        surf(ax2,x,z(indn),obs.energyDensityBox(:,indn,i1,2)');
        view(2); shading flat; box on;
        cb2 = colorbar; clim(ax2,[0 cmax(2)]); 
        cmap = colormap(ax2,'pink'); colormap(ax2,cmap(end:-1:1,:)); 
        set(ax2,'XLim',[min(x) max(x)],'YLim',[min(z(indn)) max(z(indn))], ...
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
function plotTotalEnergy_old( obs, cmax )

% constants
Nt = length(obs.t);
x = obs.boxX;
z = obs.boxZ;

% estimate cmax and make sure low values are not plotted
if isempty(cmax)
    cmax = squeeze(max(max(obs.energyDensityBox,[],1),[],2));
    cmax = cmax(ceil(obs.Nt/2),:);
end
level1  = logspace(-4,log10(cmax(1)),10);
level2  = logspace(-4,log10(cmax(2)),10);

% plot total energy - loop on time

for i1=1:Nt
    if i1==1; figure; lappend = false; else; clf; lappend = true; end
    ax1 = axes;
%    hold off; contour(ax1,x,z,obs.energyDensityBox(:,:,i1,1)',level1);
    hold off; surf(ax1,x,z,obs.energyDensityBox(:,:,i1,1)');
    view(2); shading flat;% alpha 0.5;
    colormap(ax1,'pink'); clim(ax1,[0 cmax(1)])
    set(ax1,'Position',[.15 .11 .685 .815],'colorscale','log');
    ax1.XLim = [min(x) max(x)]; ax1.YLim = [min(z) max(z)];
    ax1.PlotBoxAspectRatio = [range(x) range(z) 1];
    grid(ax1,'off')
    cb1 = colorbar(ax1,'Position',[.84 .11 .03 .815]);
    title(ax1,['Total Energy Density, time t = ' num2str(obs.t(i1)) 's'])
    if ~obs.acoustics
        title(cb1,'P energy')
        ax2 = axes;
        surf(ax2, x,z,obs.energyDensityBox(:,:,i1,2)');
%         contour(ax2,x,z,obs.energyDensityBox(:,:,i1,2)',level2);
       view(2); shading flat; alpha 0.5
        cmap = colormap(ax2,'pink'); colormap(ax2,cmap(end:-1:1,:)); 
        clim(ax2,[0 cmax(2)]);
        set(ax2,'colorscale','log');
        ax2.XLim = [min(x) max(x)]; ax2.YLim = [min(z) max(z)];
        ax2.PlotBoxAspectRatio = [range(x) range(z) 1];
        linkaxes([ax1,ax2])
        set([ax1,ax2],'Position',[.1 .11 .685 .815]);
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        cb1.Position = [.8 .11 .03 .815];
        cb2 = colorbar(ax2,'Position',[.9 .11 .03 .815]);
        title(cb2,'S energy')
    end

    % export graphics
    exportgraphics(gcf,'movieTotalEnergy.gif','Append',lappend);
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
obsEc1 = obs.energyDomainCoherent(:,1);
obsEi1 = obs.energyIncoherent(:,ind,:,1);
if ~obs.acoustics
    obsEc2 = obs.energyDomainCoherent(:,2);
    obsEi2 = obs.energyIncoherent(:,ind,:,2);
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
