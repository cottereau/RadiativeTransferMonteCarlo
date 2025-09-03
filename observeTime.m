function E = observeTime(geometry,acoustics,x,p,dir,bins,ibins,vals)

% constants
d = geometry.dimension;
frame = geometry.frame;
N1 = length(bins{1})-1;
N2 = length(bins{2})-1;
E = zeros(N1,N2,1,1+~acoustics,'uint32');

% compute radius and position angles
if d==2
    [theta,r] = cart2pol(x(:,1),x(:,2));
else
    [theta,phi,r] = cart2sph(x(:,1),x(:,2),x(:,3));
end

% compute angle between direction of propagation and position
% check normalization of dir
cospsi = dot(x,dir,2)./r;
cospsi(cospsi>1)=1;
cospsi(cospsi<-1)=-1;
% for 3D, the direction bins correspond to cos(psi)
if d==3
    psi = cospsi;
else
    psi = acos(cospsi);
end

% 3D spherical
if d==3 && strcmp(frame,'cylindrical')
    [theta,r] = cart2pol(x(:,1),x(:,2));
end

% in 3D spherical, binZ corresponds to sin(phi)
if d==3 && strcmp(frame,'spherical')
    phi = sin(phi);
end

% prepare histograms
if strcmp(frame,'cartesian')
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
elseif strcmp(frame,'cylindrical')
    if all(ibins==[1 2])
        int1 = psi;
        if d==3, int2 = x(:,3); end
        hist1 = r;
        hist2 = theta;
    elseif all(ibins==[1 3])
        int1 = theta;
        int2 = psi;
        hist1 = r;
        hist2 = x(:,3);
    elseif all(ibins==[1 4])
        int1 = theta;
        if d==3, int2 = x(:,3); end
        hist1 = r;
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
