function obs = observePaths(obs,ip,it,P)

if ip==1 && obs.check
    Np = size(obs.path,3);
    obs.path(it,:,:) = reshape([P.x(1:Np) P.y(1:Np) P.p(1:Np)]',[1 3 Np]);
end

end
