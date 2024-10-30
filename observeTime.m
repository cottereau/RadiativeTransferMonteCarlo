function E = observeTime(geometry,acoustics,x,p,dir,bins,ibins,vals)

% constants
d = geometry.dimension;
N1 = length(bins{1})-1;
N2 = length(bins{2})-1;
E = zeros(N1,N2,1,1+~acoustics,'uint32');

% compute radius and position angles
if d==2
    [theta,r] = cart2pol(x(:,1),x(:,2));
elseif d==3
    [theta,phi,r] = cart2sph(x(:,1),x(:,2),x(:,3));
end

% compute angle between direction of propagation and position
% check normalization of dir
cospsi = dot(x,dir,2)./r;
cospsi(cospsi>1)=1;
cospsi(cospsi<-1)=-1;
% for d==3, the bins correspond to cos(psi)
if d==2
    psi = acos(cospsi);
elseif d==3
    psi = cospsi;
end

% prepare histograms
if isfield(geometry,'frame') && strcmp(geometry.frame,'cartesian')
    if all(ibins==[1 2])
        int1 = psi;
        if d==3, int2 = x(:,3); end
        hist1 = x(:,1);
        hist2 = x(:,2);
    elseif all(ibins==[1 3])
        int1 = x(:,2);
        int2 = psi;
        hist1 = x(:,1);
        hist2 = x(:,3);
    elseif all(ibins==[1 4])
        int1 = x(:,2);
        if d==3, int2 = x(:,3); end
        hist1 = x(:,1);
        hist2 = psi;
    end
else
    if all(ibins==[1 2])
        int1 = psi;
        if d==3, int2 = phi; end
        hist1 = r;
        hist2 = theta;
    elseif all(ibins==[1 3])
        int1 = theta;
        int2 = psi;
        hist1 = r;
        hist2 = phi;
    elseif all(ibins==[1 4])
        int1 = theta;
        if d==3, int2 = phi; end
        hist1 = r;
        hist2 = psi;
    end
end

% accumulate energies
ind = int1>=vals{1}(1) & int1<=vals{1}(2);
if d==3 || all(ibins==[1 3])
    ind = ind & int2>=vals{2}(1) & int2<=vals{2}(2);
end
E(:,:,1,1) = histcounts2(  hist1(p&ind), hist2(p&ind), bins{1}, bins{2} );
if ~acoustics
    E(:,:,1,2) = histcounts2( hist1(~p&ind), hist2(~p&ind), bins{1}, bins{2} );
end


% % energy as a function of (r,theta) for fixed values of phi and psi
% if all(ibins==[1 2])
%     ind = psi>=vals{2}(1) & psi<=vals{2}(2);
%     if d==3
%         ind = ind & phi>=vals{1}(1) & phi<=vals{1}(2);
%     end
%     E(:,:,1,1) = histcounts2(  r(p&ind), theta(p&ind), bins{1}, bins{2} );
%     if ~acoustics
%         E(:,:,1,2) = histcounts2( r(~p&ind), theta(~p&ind), bins{1}, bins{2} );
%     end
%     % energy as a function of (r,phi) for fixed values of theta and psi
% elseif all(ibins==[1 3])
%     ind = theta>=vals{1}(1) & theta<=vals{1}(2) ...
%         & psi>=vals{2}(1) & psi<=vals{2}(2);
%     E(:,:,1,1) = histcounts2(  r(p&ind), phi(p&ind), bins{1}, bins{2} );
%     if ~acoustics
%         E(:,:,1,2) = histcounts2( r(~p&ind), phi(~p&ind), bins{1}, bins{2} );
%     end
%     % energy as a function of (r,psi) for fixed values of theta and phi
% elseif all(ibins==[1 4])
%     ind = theta>=vals{1}(1) & theta<=vals{1}(2);
%     if d==3
%         ind = ind & phi>=vals{2}(1) & phi<=vals{2}(2);
%     end
%     E(:,:,1,1) = histcounts2(  r(p&ind), psi(p&ind), bins{1}, bins{2} );
%     if ~acoustics
%         E(:,:,1,2) = histcounts2( r(~p&ind), psi(~p&ind), bins{1}, bins{2} );
%     end
% end
