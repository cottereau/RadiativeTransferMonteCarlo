function scatterDirections(obs,station)

% dimensions
Nt = length(obs.t);

figure;
for i1 = 1:length(station)
    v = VideoWriter(['scatterPlot' num2str(station(i1)) '.avi']);
    open(v);
    for i2=1:Nt
        pol = obs.energy(:,station(i1),i2);
        polarplot(obs.psi,pol);
        title(['time T = ' num2str(obs.t(i2)) 's'])
        M(i2) = getframe;
        writeVideo(v,M(i2));
    end
    close(v);
end
