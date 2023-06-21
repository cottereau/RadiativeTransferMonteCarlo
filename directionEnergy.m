function [psi2pi,Ec,Ei] = directionEnergy(obs,sensors)
% construct directional energy at sensors

% constant
Ns = size(sensors,1);
psi2pi = [2*pi-obs.psi(end:-1:2) obs.psi];

% initialization
Ec = zeros(2*obs.Npsi-1,obs.Nt,Ns);
Ei = zeros(2*obs.Npsi-1,obs.Nt,Ns);
Ecref = obs.energy(:,:,:,1);
Eiref = obs.energy(:,:,:,2);

% loop on sources
for i2 = 1:obs.nSources

    % loop on sensors
    for i1 = 1:Ns

        % distance and angle from source to sensor
        X = sensors(i1,:) - obs.positionSources(i2,:);
        r = sqrt(sum(X.^2));
        rxz = sqrt(sum(X([1 3]).^2));
        theta = acos(X(1)/rxz);
        if X(3)<0; theta = 2*pi-theta; end
        psiRot = mod( psi2pi - theta, 2*pi );

        % estimate coherent directional energy at sensor and rotate
        Es = squeeze( interp1( obs.r', Ecref, r ));
        Es = [ Es(end:-1:2,:); Es];
        Ec(:,:,i1) = Ec(:,:,i1) + interp1( psi2pi', Es, psiRot );

        % estimate incoherent directional energy at sensor and rotate
        Es = squeeze( interp1( obs.r', Eiref, r ));
        Es = [ Es(end:-1:2,:); Es];
        Ei(:,:,i1) = Ei(:,:,i1) + interp1( psi2pi', Es, psiRot );

    % end of loop on sources
    end

% end of loop on sensors
end
