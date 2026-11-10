function plotEnergies( type, obs, material, lambda, cmax, rmax )
% plotEnergies  Plot or export energy diagnostics from radiativeTransfer.
%
% For total-energy movies, integrate observations over all propagation
% directions, so that observation.directions = [0 pi] and obs.Npsi == 1.
% For elastic simulations, the total-energy movie displays three linked
% panels containing P, S, and P+S energy.
%
% Optional fields of TYPE for total-energy movies:
%   energyScale        'linear' (default) or 'log10'.
%   clim               Fixed two-element color limits shared by all panels,
%                      or a 3-by-2 array for elastic P, S, and total panels.
%                      By default, linear limits use the maximum energy at
%                      the middle observation time for each panel.
%   logFloor           Positive floor used by log10. Default: 1e-12.
%   colormap           MATLAB colormap name. Default: 'turbo'.
%   filename           Output GIF filename. Default: movieTotalEnergy.gif.
%   delayTime          Delay between GIF frames. Default: 0.05 seconds.
%   visible            Figure visibility, 'on' (default) or 'off'.
%   axisMode           'normal' (default), 'equal', or 'tight'.
%   mirrorCylindrical  Mirror cylindrical r-z plots across r = 0. The
%                      default is true for cylindrical r-z plots.
%   xlim, zlim         Optional fixed horizontal and vertical plot limits.
%   boundaries         Optional geometry.bnd array to overlay.
%
% Example:
%   plotting = struct('movieTotalEnergy',true, ...
%                     'energyScale','log10', ...
%                     'clim',[-8 -2], ...
%                     'filename','movieTotalEnergy.gif', ...
%                     'visible','off');
%   plotEnergies(plotting, obs, material, source.lambda);

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
        plotTotalEnergy( obs, cmax, type )
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
function plotTotalEnergy(obs,cmax,type)

% Identify the two spatial coordinates stored in energyDensity.
if obs.Nx == 1
    horizontalValues = obs.y;
    verticalValues = obs.z;
    horizontalEdges = obs.binY;
    verticalEdges = obs.binZ;
    horizontalCoordinate = 2;
    verticalCoordinate = 3;
elseif obs.Ny == 1
    horizontalValues = obs.x;
    verticalValues = obs.z;
    horizontalEdges = obs.binX;
    verticalEdges = obs.binZ;
    horizontalCoordinate = 1;
    verticalCoordinate = 3;
else
    horizontalValues = obs.x;
    verticalValues = obs.y;
    horizontalEdges = obs.binX;
    verticalEdges = obs.binY;
    horizontalCoordinate = 1;
    verticalCoordinate = 2;
end

horizontalPlot = plottingCoordinates(horizontalValues,horizontalEdges, ...
    size(obs.energyDensity,1));
verticalPlot = plottingCoordinates(verticalValues,verticalEdges, ...
    size(obs.energyDensity,2));

if isfield(obs,'frame')
    frame = obs.frame;
else
    frame = 'cartesian';
    if isfield(type,'boundaries') && ~isempty(type.boundaries) && ...
            any([type.boundaries.dir] == 4)
        frame = 'cylindrical';
    end
end

horizontalLabel = coordinateLabel(frame,horizontalCoordinate);
verticalLabel = coordinateLabel(frame,verticalCoordinate);
isCylindricalRadius = strcmp(frame,'cylindrical') && ...
    horizontalCoordinate == 1;
isCylindricalRZ = isCylindricalRadius && verticalCoordinate == 3;

% Read plotting options and apply defaults.
if ~isfield(type,'energyScale') || isempty(type.energyScale)
    type.energyScale = 'linear';
end
if ~isfield(type,'colormap') || isempty(type.colormap)
    type.colormap = 'turbo';
end
if ~isfield(type,'filename') || isempty(type.filename)
    type.filename = 'movieTotalEnergy.gif';
end
if ~isfield(type,'delayTime') || isempty(type.delayTime)
    type.delayTime = 0.05;
end
if ~isfield(type,'logFloor') || isempty(type.logFloor)
    type.logFloor = 1e-12;
end
if ~isfield(type,'visible') || isempty(type.visible)
    type.visible = 'on';
end
if ~isfield(type,'axisMode') || isempty(type.axisMode)
    type.axisMode = 'normal';
end
if ~isfield(type,'mirrorCylindrical') || isempty(type.mirrorCylindrical)
    type.mirrorCylindrical = isCylindricalRZ;
end

if ~isscalar(type.mirrorCylindrical) || ...
        ~(islogical(type.mirrorCylindrical) || isnumeric(type.mirrorCylindrical))
    error('type.mirrorCylindrical must be a scalar logical value.');
end
mirrorCylindrical = logical(type.mirrorCylindrical);
if mirrorCylindrical && ~isCylindricalRZ
    error('Cylindrical mirroring is only available for cylindrical r-z plots.');
end
if ~isscalar(type.delayTime) || ~isfinite(type.delayTime) || type.delayTime < 0
    error('type.delayTime must be a finite nonnegative scalar.');
end
if ~isscalar(type.logFloor) || ~isfinite(type.logFloor) || type.logFloor <= 0
    error('type.logFloor must be a finite positive scalar.');
end

% Assemble the acoustic or elastic energy panels.
energyDensity = obs.energyDensity;
if obs.acoustics
    energyPanels = energyDensity(:,:,:,1);
    panelTitles = {'Energy density'};
else
    pressureEnergy = energyDensity(:,:,:,1);
    shearEnergy = energyDensity(:,:,:,2);
    energyPanels = cat(4,pressureEnergy,shearEnergy, ...
        pressureEnergy + shearEnergy);
    panelTitles = {'P energy density','S energy density', ...
        'Total energy density'};
end
numberPanels = numel(panelTitles);

if mirrorCylindrical
    horizontalPlot = [-horizontalPlot(end:-1:1) horizontalPlot];
    energyPanels = cat(1,energyPanels(end:-1:1,:,:,:),energyPanels);
    horizontalLabel = 'x';
end

% Determine color scaling.
switch lower(char(type.energyScale))
    case 'linear'
        plottedEnergy = energyPanels;
        if isfield(type,'clim') && ~isempty(type.clim)
            colorLimits = validateColorLimits(type.clim,numberPanels);
        elseif ~isempty(cmax)
            if ~isnumeric(cmax) || any(~isfinite(cmax(:))) || any(cmax(:) <= 0)
                error('cmax must contain finite positive values.');
            end
            if numberPanels == 1
                maximumEnergy = max(cmax(:));
            elseif isscalar(cmax)
                maximumEnergy = repmat(cmax,1,numberPanels);
            elseif numel(cmax) == 2
                maximumEnergy = [cmax(:)' sum(cmax(:))];
            elseif numel(cmax) == numberPanels
                maximumEnergy = cmax(:)';
            else
                error('For elastic movies, cmax must contain one, two, or three values.');
            end
            colorLimits = [zeros(numberPanels,1) maximumEnergy(:)];
        else
            % Use the middle observation time as a representative fixed
            % linear scale, following the historical plotting behavior.
            referenceTimeIndex = ceil(obs.Nt/2);
            colorLimits = zeros(numberPanels,2);
            for panelIndex = 1:numberPanels
                referenceEnergy = energyPanels(:,:,referenceTimeIndex,panelIndex);
                finiteReferenceEnergy = referenceEnergy(isfinite(referenceEnergy));

                if isempty(finiteReferenceEnergy)
                    maximumEnergy = 0;
                else
                    maximumEnergy = max(finiteReferenceEnergy);
                end

                % Fall back to the global panel maximum if the reference
                % frame contains no finite positive energy.
                if maximumEnergy <= 0
                    panelEnergy = energyPanels(:,:,:,panelIndex);
                    finitePanelEnergy = panelEnergy(isfinite(panelEnergy));
                    if isempty(finitePanelEnergy)
                        maximumEnergy = 1;
                    else
                        maximumEnergy = max(finitePanelEnergy);
                        if maximumEnergy <= 0
                            maximumEnergy = 1;
                        end
                    end
                end
                colorLimits(panelIndex,:) = [0 maximumEnergy];
            end
        end

    case 'log10'
        plottedEnergy = log10(max(energyPanels,type.logFloor));
        if isfield(type,'clim') && ~isempty(type.clim)
            colorLimits = validateColorLimits(type.clim,numberPanels);
        else
            colorLimits = repmat([-8 -2],numberPanels,1);
        end

    otherwise
        error('Unknown energyScale "%s". Use "linear" or "log10".', ...
            char(type.energyScale));
end

% Determine fixed spatial limits.
if isfield(type,'xlim') && ~isempty(type.xlim)
    horizontalLimits = validateAxisLimits(type.xlim,'type.xlim');
    if mirrorCylindrical && horizontalLimits(1) >= 0
        horizontalLimits = [-horizontalLimits(2) horizontalLimits(2)];
    end
else
    horizontalLimits = defaultAxisLimits(horizontalPlot,horizontalEdges);
    if mirrorCylindrical
        horizontalLimits = [-max(abs(horizontalLimits)) max(abs(horizontalLimits))];
    end
end

if isfield(type,'zlim') && ~isempty(type.zlim)
    verticalLimits = validateAxisLimits(type.zlim,'type.zlim');
else
    verticalLimits = defaultAxisLimits(verticalPlot,verticalEdges);
end

% Replace an earlier movie instead of appending to it.
filename = char(type.filename);
if exist(filename,'file')
    delete(filename);
end

fig = figure('Color','w','Visible',type.visible);
layout = tiledlayout(fig,numberPanels,1, ...
    'TileSpacing','compact','Padding','compact');
axesHandles = gobjects(numberPanels,1);
for panelIndex = 1:numberPanels
    axesHandles(panelIndex) = nexttile(layout,panelIndex);
end
linkaxes(axesHandles,'xy');

for timeIndex = 1:length(obs.t)
    for panelIndex = 1:numberPanels
        ax = axesHandles(panelIndex);
        cla(ax);
        imagesc(ax,horizontalPlot,verticalPlot, ...
            plottedEnergy(:,:,timeIndex,panelIndex)');
        set(ax,'YDir','normal');
        colormap(ax,type.colormap);
        clim(ax,colorLimits(panelIndex,:));
        colorbar(ax);

        switch lower(char(type.axisMode))
            case 'equal'
                axis(ax,'equal');
            case 'tight'
                axis(ax,'tight');
            case 'normal'
                axis(ax,'normal');
            otherwise
                error('Unknown axisMode "%s". Use "normal", "equal", or "tight".', ...
                    char(type.axisMode));
        end

        xlim(ax,horizontalLimits);
        ylim(ax,verticalLimits);
        box(ax,'on');
        ylabel(ax,verticalLabel);
        title(ax,sprintf('%s, t = %.3g s', ...
            panelTitles{panelIndex},obs.t(timeIndex)));

        if panelIndex == numberPanels
            xlabel(ax,horizontalLabel);
        else
            set(ax,'XTickLabel',[]);
        end

        if isfield(type,'boundaries') && ~isempty(type.boundaries)
            hold(ax,'on');
            plotBoundaries(ax,type.boundaries,frame,horizontalCoordinate, ...
                verticalCoordinate,isCylindricalRadius,mirrorCylindrical);
            hold(ax,'off');
        end
    end

    drawnow;
    movieFrame = getframe(fig);
    [imageData,colorMap] = rgb2ind(frame2im(movieFrame),256);
    if timeIndex == 1
        imwrite(imageData,colorMap,filename,'gif','LoopCount',Inf, ...
            'DelayTime',type.delayTime);
    else
        imwrite(imageData,colorMap,filename,'gif','WriteMode','append', ...
            'DelayTime',type.delayTime);
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coordinates = plottingCoordinates(values,edges,numberValues)
if numel(values) == numberValues
    coordinates = values;
elseif numel(edges) == numberValues + 1
    coordinates = (edges(1:end-1) + edges(2:end))/2;
else
    error('Observation coordinates do not match the energy-density dimensions.');
end
coordinates = coordinates(:)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label = coordinateLabel(frame,coordinate)
labels = {'x','y','z'};
if strcmp(frame,'cylindrical')
    labels = {'r','azimuth','z'};
elseif strcmp(frame,'spherical')
    labels = {'r','azimuth','elevation'};
end
label = labels{coordinate};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limits = defaultAxisLimits(values,edges)
if numel(edges) >= 2 && all(isfinite(edges(:)))
    limits = [min(edges(:)) max(edges(:))];
else
    limits = [min(values(:)) max(values(:))];
end
if limits(2) <= limits(1)
    padding = max(0.5,0.5*abs(limits(1)));
    limits = limits(1) + [-padding padding];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limits = validateAxisLimits(limits,name)
if ~isnumeric(limits) || numel(limits) ~= 2 || ...
        any(~isfinite(limits)) || limits(2) <= limits(1)
    error('%s must contain two finite increasing values.',name);
end
limits = limits(:)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function limits = validateColorLimits(limits,numberPanels)
if isnumeric(limits) && isvector(limits) && numel(limits) == 2
    limits = validateAxisLimits(limits,'type.clim');
    limits = repmat(limits,numberPanels,1);
elseif isnumeric(limits) && isequal(size(limits),[numberPanels 2])
    for panelIndex = 1:numberPanels
        limits(panelIndex,:) = validateAxisLimits( ...
            limits(panelIndex,:),'each row of type.clim');
    end
else
    error('type.clim must be a two-element vector or a %d-by-2 array.', ...
        numberPanels);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotBoundaries(ax,boundaries,frame,horizontalCoordinate, ...
        verticalCoordinate,isCylindricalRadius,mirrorCylindrical)
for boundaryIndex = 1:numel(boundaries)
    boundary = boundaries(boundaryIndex);
    if strcmp(frame,'cartesian')
        if boundary.dir == horizontalCoordinate
            xline(ax,boundary.val,'w--','LineWidth',1.5);
        elseif boundary.dir == verticalCoordinate
            yline(ax,boundary.val,'w--','LineWidth',1.5);
        end
    elseif strcmp(frame,'cylindrical')
        if boundary.dir == 3 && verticalCoordinate == 3
            yline(ax,boundary.val,'w--','LineWidth',1.5);
        elseif boundary.dir == 4 && isCylindricalRadius
            xline(ax,boundary.val,'w--','LineWidth',1.5);
            if mirrorCylindrical
                xline(ax,-boundary.val,'w--','LineWidth',1.5);
            end
        end
    end
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
        warning('%s',er)
    end
end

end
