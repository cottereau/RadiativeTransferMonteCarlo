function E = observeTime(d,acoustics,x,p,dir,bins,ibins,vals)

% constants
Nt = size(x,3);
N1 = length(bins{1})-1;
N2 = length(bins{2})-1;
E = zeros(N1,N2,Nt,1+~acoustics,'uint32');

% loop on time steps
for it = 1:Nt

    % compute radius and position angles
    if d==2
        [theta,r] = cart2pol(x(:,1,it),x(:,2,it));
    elseif d==3
        [theta,phi,r] = cart2sph(x(:,1,it),x(:,2,it),x(:,3,it));
    end

    % compute angle between direction of propagation and position
    % check normalization of dir
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
    % energy as a function of r x theta for fixed values of phi and psi
    if all(ibins==[1 2])
        ind = psi>=vals{2}(1) & psi<=vals{2}(2);
        if d==3
            ind = ind & phi>=vals{1}(1) & phi<=vals{1}(2);
        end
        E(:,:,it,1) = histcounts2(  r(p(:,it)&ind), theta(p(:,it)&ind), bins{1}, bins{2} );
        if ~acoustics
            E(:,:,it,2) = histcounts2( r(~p(:,it)&ind), theta(~p(:,it)&ind), bins{1}, bins{2} );
        end
    % energy as a function of r x phi for fixed values of theta and psi
    elseif all(ibins==[1 3])
        ind = theta>=vals{1}(1) & theta<=vals{1}(2) ...
                                       & psi>=vals{2}(1) & psi<=vals{2}(2);
        E(:,:,it,1) = histcounts2(  r(p(:,it)&ind), phi(p(:,it)&ind), bins{1}, bins{2} );
        if ~acoustics
            E(:,:,it,2) = histcounts2( r(~p(:,it)&ind), phi(~p(:,it)&ind), bins{1}, bins{2} );
        end
    % energy as a function of r x psi for fixed values of theta and phi
    elseif all(ibins==[1 4])
        ind = theta>=vals{1}(1) & theta<=vals{1}(2);
        if d==3
            ind = ind & phi>=vals{2}(1) & phi<=vals{2}(2);
        end
        E(:,:,it,1) = histcounts2(  r(p(:,it)&ind), psi(p(:,it)&ind), bins{1}, bins{2} );
        if ~acoustics
            E(:,:,it,2) = histcounts2( r(~p(:,it)&ind), psi(~p(:,it)&ind), bins{1}, bins{2} );
        end
    end
end
