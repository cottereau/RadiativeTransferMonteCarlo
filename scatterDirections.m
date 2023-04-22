function scatterDirections(obs,rmax)

% dimensions
Nt = length(obs.t);
Ns = size(obs.sensors,1);

figure;
lappend = false;
for i2 = 1:Nt
    for i1 = 1:Ns
        subplot(1,Ns,i1)
        polarplot(obs.psi2pi,obs.energyDirectional(:,i2,i1));
        rlim([0 rmax])
    end
    title(['time T = ' num2str(obs.t(i2)) 's'])
    exportgraphics(gcf,'testScatter.gif','Append',lappend);
    lappend = true;
end
