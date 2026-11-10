function [measure,integrated] = observationCellMeasure(obs,axes)
% Full spatial-cell measure, with an angular factor only if resolved.
% A single infinite-width spatial bin is an integral/marginal, not a
% finite-volume density. Record that axis and do not divide by infinity.
widths = {obs.dx,obs.dy,obs.dz,obs.dpsi};
measure = ones(numel(widths{axes(1)}),numel(widths{axes(2)}));
integrated = [];
used = 1:obs.d;
if any(axes==4), used(end+1) = 4; end
for k = used
    width = widths{k};
    if any(isinf(width))
        if numel(width)~=1
            error('RTE:InfiniteResolvedCell','Infinite spatial bins can only be used as a single integrated interval.');
        end
        integrated(end+1) = k; %#ok<AGROW>
        width = 1;
    end
    if any(~isfinite(width)) || any(width<=0)
        error('RTE:CellMeasure','Cell measures must be positive and finite.');
    end
    if k==axes(1)
        measure = measure.*width(:);
    elseif k==axes(2)
        measure = measure.*width(:)';
    else
        measure = measure*width;
    end
end
end
