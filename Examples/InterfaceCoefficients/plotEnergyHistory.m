function plotEnergyHistory(ax, out, expected, keep, titleStr)
%PLOTENERGYHISTORY  Energy vs time for one interface case, with the level
%   each channel is supposed to settle at drawn over it.
%
%   plotEnergyHistory(ax, out, expected, keep, titleStr) plots the channels
%   listed in keep (indices into out.E's columns) as solid curves against
%   t / tHitCenter, the total as a dashed black line, and the theoretical
%   final level of each kept channel as a dotted horizontal line in the same
%   colour. The gap between a curve's plateau and its dotted line is the
%   error the summary table reports numerically.
%
%   out       as returned by interfaceEnergyRTE
%   expected  1 x 7 theoretical final energies (NaN where not applicable)
%   keep      indices of the channels worth drawing for this case
%
%   See also INTERFACEENERGYRTE, PLOTCHANNELS.

cols = lines(7);
tn = out.t / out.tHitCenter;

hold(ax,'on');
h = gobjects(0); names = {};
for k = keep
    plot(ax, tn, out.E(:,k), '-', 'Color', cols(k,:), 'LineWidth', 1.4);
    h(end+1) = plot(ax, NaN, NaN, '-', 'Color', cols(k,:), 'LineWidth', 1.4); %#ok<AGROW>
    names{end+1} = strrep(out.labels{k},'_','\_'); %#ok<AGROW>
    if isfinite(expected(k))
        yline(ax, expected(k), ':', 'Color', cols(k,:), 'LineWidth', 1.2);
    end
end
plot(ax, tn, out.E(:,7), 'k--', 'LineWidth', 1.1);
h(end+1) = plot(ax, NaN, NaN, 'k--', 'LineWidth', 1.1); names{end+1} = 'total';
xline(ax, 1, 'k:', 'LineWidth', 0.8);

xlabel(ax,'t / t_{hit,center}');
ylabel(ax,'energy fraction');
title(ax, titleStr);
ylim(ax,[-0.05 1.05]);
grid(ax,'on'); box(ax,'on');
legend(ax, h, names, 'Location','best', 'FontSize', 7);
end
