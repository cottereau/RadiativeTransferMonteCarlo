function plotChannels(figTitle, theta_grid, theoryCurves, labels, testAngles, mcValues)
%PLOTCHANNELS  Theory curve (line) vs Monte Carlo (markers), per channel.
%
%   plotChannels(figTitle, theta_grid, theoryCurves, labels, testAngles, mcValues)
%   mirrors BoundaryReflectionEnergyTests.m's validation-plot idea -- code
%   output against the theoretical reference -- adapted to this folder's
%   angle-sweep format: BoundaryReflectionEnergyTests.m plots energy vs
%   TIME for one fixed case, so its theoretical reference is a constant
%   and draws as a flat horizontal line; these examples instead sweep
%   incidence ANGLE for several channels at once, so the natural analogue
%   of "the line the code's output should sit on" is the full theoretical
%   R/T curve over angle (theoryCurves, evaluated on theta_grid, e.g. the
%   fine j1_deg grid these examples already build for interp1), with the
%   Monte Carlo fractions at just the tested angles (testAngles, mcValues)
%   marked on top of it. A converged simulation's markers sit ON the
%   curve, the same visual check a flat theory line gives.
%
%   theoryCurves, mcValues: cell arrays, one entry per channel, matching
%   labels (one legend entry per channel, theory and MC sharing colour).

figure('Name', figTitle);
hold on
co = get(gca,'ColorOrder');
for c = 1:numel(labels)
    col = co(mod(c-1,size(co,1))+1,:);
    plot(theta_grid, theoryCurves{c}, '-', 'Color', col, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%s (theory)', labels{c}));
    plot(testAngles{c}, mcValues{c}, 'o', 'Color', col, 'MarkerFaceColor', col, ...
        'MarkerSize', 6, 'DisplayName', sprintf('%s (MC)', labels{c}));
end
hold off
grid on
xlabel('incidence angle (deg)');
ylabel('energy fraction');
ylim([-0.05 1.05]);
legend('Location','best');
title(figTitle, 'Interpreter', 'none');
end
