function scatterDirections(obs,station)
figure;
Nt = length(obs.t);
theta = [-obs.d(end:-1:1) obs.d];
for i1 = 1:length(station)
    v = VideoWriter(['scatterPlot' num2str(station(i1)) '.avi']);
    open(v);
    for i2=1:Nt
        pol = [obs.energy(end:-1:1,i2,station(i1)); ...
                                     obs.energy(:,i2,station(i1)) ];
        polarplot(theta,pol);
        title(['time T = ' num2str(obs.t(i2)) 's'])
        M(i2) = getframe;
        writeVideo(v,M(i2));
    end
    close(v);
end
