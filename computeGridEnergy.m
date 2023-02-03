function obs = computeGridEnergy(obs,acoustics)
Np = length(obs.x);
Nt = length(obs.t);
L = max(obs.x)/sqrt(2);
xg = linspace(0,L,min(4*Np,100));
xg = [-xg(end:-1:2) xg];
Ng = length(xg);
[X,Y] = ndgrid(xg,xg);
rg = sqrt(X.^2+Y.^2);
obs.grid = xg;
obs.gridEnergy = zeros(Ng,Ng,Nt);
for i1 = 1:Nt
    obs.gridEnergy(:,:,i1) = interp1(obs.x,obs.energyDensity(:,i1),rg);
end
if ~acoustics
    obs.gridEnergyP = zeros(length(xg),length(xg),Nt);
    for i1 = 1:Nt
        obs.gridEnergyP(:,:,i1) = interp1(obs.x,obs.sensorEnergyP(:,i1),rg);
    end
    obs.gridEnergyS = zeros(length(xg),length(xg),Nt);
    for i1 = 1:Nt
        obs.gridEnergyS(:,:,i1) = interp1(obs.x,obs.sensorEnergyS(:,i1),rg);
    end
end
end