function fig = plotSourceBallisticComparison(result)
%PLOTSOURCEBALLISTICCOMPARISON RTE and analytic cell energies at three times.
% Each x tick is a radial/elevation cell; values are fractions of modal J.
fig=figure('Name',result.caseID);
modes={'P','S'};
for beta=1:2
    channel=result.(modes{beta});
    for it=1:3
        subplot(2,3,(beta-1)*3+it);
        if isempty(channel)
            axis off; title('No S radiation'); continue
        end
        measured=channel.measuredCellEnergyJ(:,:,it)/channel.energyJ;
        analytic=channel.analyticCellEnergyJ(:,:,it)/channel.energyJ;
        finitePulse=channel.finitePulseFarFieldCellEnergyJ(:,:,it)/channel.energyJ;
        bar([analytic(:) measured(:) finitePulse(:)]);
        xlabel('Radial/elevation cell'); ylabel('Fraction of modal energy');
        title(sprintf('%s, t = %.4g s',modes{beta},result.observation.time(it)));
        if beta==1 && it==1
            legend('Event analytic','RTE','Finite pulse far field');
        end
    end
end
end
