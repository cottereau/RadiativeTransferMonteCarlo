function scatterDirections(obs,station)

% dimensions
Nt = length(obs.t);
Nphi = obs.Nphi;

figure;
for i1 = 1:length(station)
    v = VideoWriter(['scatterPlot' num2str(station(i1)) '.avi']);
    open(v);
    for i2=1:Nt
        pol = obs.energy(station(i1),:,i2);
        ind = 1:floor(Nphi/2);
        pol(ind) = (pol(ind)+pol(Nphi+1-ind))/2;
        pol(Nphi+1-ind) = pol(ind);
        polarplot([obs.phi obs.phi(1)],[pol pol(1)]);
        title(['time T = ' num2str(obs.t(i2)) 's'])
        M(i2) = getframe;
        writeVideo(v,M(i2));
    end
    close(v);
end
