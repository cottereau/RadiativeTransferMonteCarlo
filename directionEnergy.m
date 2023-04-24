function [psi2pi,E] = directionEnergy(obs,sensors)
% construct directional energy at sensors

% constant
Ns = size(sensors,1);
psi2pi = [0 obs.psi 2*pi-obs.psi(end:-1:1) 2*pi];

% initialization
E = zeros(2*(obs.Npsi+1),obs.Nt,Ns);

% loop on sensors
for i1 = 1:Ns

    % loop on sources
    for i2 = 1:obs.nSources

        % distance and angle from source to sensor
        X = sensors(i1,:) - obs.positionSources(i2,:);
        r = sqrt(sum(X.^2));
        rxz = sqrt(sum(X([1 3]).^2));
        theta = acos(X(1)/rxz);
        if X(3)<0; theta = 2*pi-theta; end

        % estimate directional energy at sensor
        Es = squeeze( interp1( obs.r', permute( obs.energy,[2 1 3]), r ));
        Es = [ Es(1,:); Es; Es(end,:); Es(end:-1:1,:) ];

        % rotate the directional energy
        psiRot = mod( psi2pi - theta, 2*pi );
        E(:,:,i1) = E(:,:,i1) + interp1( psi2pi', Es, psiRot );

    % end of loop on sources
    end

% end of loop on sensors
end
