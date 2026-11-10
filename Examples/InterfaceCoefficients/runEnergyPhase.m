function allPass = runEnergyPhase(name, matArray, modeIn, thetas, expectedFun, keep, N)
%RUNENERGYPHASE  RTE energy phase shared by the InterfaceCoefficients examples.
%
%   allPass = runEnergyPhase(name, matArray, modeIn, thetas, expectedFun,
%   keep, N) runs the transport once per incidence angle with
%   interfaceEnergyRTE, prints for every channel the level the energy SHOULD
%   settle at and the level it ACTUALLY reaches, and draws one energy-vs-time
%   tile per angle with those levels overlaid.
%
%   The rest of this folder checks the interface operator in isolation, by
%   calling MaterialInterface* once on a synthetic particle set. This phase
%   is the end-to-end counterpart: particles are launched away from the
%   interface, propagate, cross it, and are measured from where they land.
%   A discrepancy that only the transport can produce -- a lost particle, a
%   mis-set matIdx, a polarization that does not survive the crossing --
%   shows up here and nowhere else.
%
%   expectedFun(theta) returns a 1x7 vector of theoretical final energies in
%   interfaceEnergyRTE's channel order, NaN where the channel cannot exist.
%   keep lists the channel indices worth printing and plotting.
%
%   See also INTERFACEENERGYRTE, PLOTENERGYHISTORY, CHECKCHANNEL.

if nargin < 7 || isempty(N), N = 50000; end
allPass = true;
tol = 5/sqrt(N);

fig = figure('Name', sprintf('%s: RTE energy histories', name), 'Color','w', ...
             'Position',[100 60 1250 820]);
nCols = 2;
tl = tiledlayout(fig, ceil(numel(thetas)/nCols), nCols, ...
                 'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s -- RTE energy vs time (dotted = level it should reach)', ...
      strrep(name,'_','\_')));

for theta = thetas
    out = interfaceEnergyRTE(matArray, modeIn, theta, N, struct('nTime',120));
    exp7 = expectedFun(theta);

    fprintf('  [RTE ] theta = %g deg   (measured total = %.4f)\n', theta, out.final(7));
    for k = keep
        allPass = checkChannel(sprintf('%s', out.labels{k}), ...
                               exp7(k), out.final(k), tol) && allPass;
    end

    plotEnergyHistory(nexttile(tl), out, exp7, keep, ...
        sprintf('%s incidence, \\theta = %g^\\circ', modeIn, theta));
end
end
