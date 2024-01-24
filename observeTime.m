function Ei = observeTime(d,acoustics,x,p,dir,coherent,binPsi,binR)

% constants
Nt = size(x,3);
Npsi = length(binPsi)-1;
Nr = length(binR)-1;
Ei = zeros(Nr,Npsi,Nt,1+~acoustics,'uint16');

% loop on time steps
for it = 1:Nt

    % compute radius and angle between direction of propagation and position
    % check normalization of dir
    r = vecnorm(x(:,:,it),2,2);
    cospsi = dot(x(:,:,it),dir(:,:,it),2)./r;
    cospsi(cospsi>1)=1;
    cospsi(cospsi<-1)=-1;

    % for d==3, the bins correspond to cos(psi)
    if d==2
        psi = acos(cospsi);
    elseif d==3
        psi = cospsi;
    end

    % accumulate energies
    Ei(:,:,it,1) = histcounts2( r(~coherent(:,it)&p(:,it)), ...
                              psi(~coherent(:,it)&p(:,it)), binR, binPsi );
    if ~acoustics
        Ei(:,:,it,2) = histcounts2( r(~coherent(:,it)&~p(:,it)), ...
                              psi(~coherent(:,it)&~p(:,it)), binR, binPsi );
    end

end
