function scatterDirections(obs,sensors)

% dimensions
Nt = length( obs.t );
Ns = size( sensors, 1 );
rmax = 1e-3;

% compute directional energy
[psi2pi,E] = directionEnergy( obs, sensors );

titl = tiledlayout(1,Ns);
lappend = false;
for i2 = 1:Nt
    for i1 = 1:Ns
        nexttile(i1)
        polarplot(psi2pi,E(:,i2,i1)); rlim([0 rmax])
        title(['sensor at [' num2str(sensors(i1,1)) ',' ...
                    num2str(sensors(i1,2)) ',' num2str(sensors(i1,3)) ']'])
    end
    title(titl,['time T = ' num2str(obs.t(i2)) 's'])
    exportgraphics(gcf,'testScatter.gif','Append',lappend);
    lappend = true;
end
