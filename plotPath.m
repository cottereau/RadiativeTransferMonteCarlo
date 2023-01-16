function plotPath(obs,source)
Np = size(obs.path,3);

% plot paths
figure;hold on
for i1=1:Np
    x = obs.path(:,1,i1);
    y = obs.path(:,2,i1);
    plot(x,y);
end

% plot source
scatter(0,0,50,'r','filled')
plotCircle(source.lambda,'r')

% plot observation circles
for i1 = 1:length(obs.rp)
    plotCircle(obs.rp(i1),'k')
end

% plot energy

function plotCircle(r,c)
theta = linspace(0,2*pi,1000);
x = r*cos(theta);
y = r*sin(theta);
 plot(x,y,[c '--'])