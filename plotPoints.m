function plotPoints(obs)
Np = length(obs.rp);
figure;plot(obs.t,obs.sensorEnergy,'linewidth',1.5)
xlabel('time [s]','FontSize',15);
ylabel('Energy','FontSize',15);
legend([repmat('x = ',[Np 1]) num2str(obs.rp)],'FontSize',15)
box on; grid on;



