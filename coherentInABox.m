function Ec = coherentInABox(Ecr,x,y,z,x0,y0,z0,t,d,lambda,c)
% Evaluate the coherent energy at positions (x,z) and time t, for a
% collection of sources at position (x0,z0) firing all at the same time
% and losing coherence with distance
%
% Ecr        amplitude of the coherent energy with distance from source
% (x,y,z)    positions at which energy is evaluated (y should be a scalar)
% (x0,y0,z0) positions of the sources
% t          times at which energy should be evaluated
% d          dimension of the problem (if d==2, y and y0 are not used)
% lambda     thickness 

% formatting
if numel(y)>1
    error('the points should be in plane (x,z) and y should be a scalar')
end
if size(x0,1)>1, x0 = x0'; end
if size(y0,1)>1, y0 = y0'; end
if size(z0,1)>1, z0 = z0'; end
if size(x,2)>1, x = x'; end
if size(y,2)>1, y = y'; end
if size(z,2)>1, z = z'; end

% constants
minR = eps;
Nt = length(t);
ns = size(x0,2);
Nx = length(x);

% prepare amplifications of energy
xs = x0-x0';
ys = y0-y0';
zs = z0-z0';
theta = cart2pol(xs,zs);
eta2 = xs.*cos(theta)+zs.*sin(theta);
phi = cart2pol(eta2,ys);
xc = (x0+x0')/2;
yc = (y0+y0')/2;
zc = (z0+z0')/2;

% preparing amplifications of coherent energy for slabs and boxes
ampC = zeros(Nx,ns);
for i1 = 1:ns
    xloc = x-xc(:,i1)';
    yloc = y-yc(:,i1)';
    zloc = z-zc(:,i1)';
    eta =  xloc.*cos(theta(:,i1)')+zloc.*sin(theta(:,i1)');
    eta =  eta.*cos(phi(:,i1)')+yloc.*sin(phi(:,i1)');
    tmpamp = 1 + exp(-((eta/lambda).^2));
    ampC(:,i1) = prod(tmpamp(:,setdiff(1:ns,i1)),2);
end

% initialization
Ec = zeros(numel(x),Nt);
if d==2
    y=0; 
    y0 = zeros(size(x0));
end

% loop on sources
for i1 = 1:length(x0)
    r = sqrt((x(:)-x0(i1)).^2+(y-y0(i1)).^2+(z(:)-z0(i1)).^2);
    r(r<minR)=minR;
    Ec = Ec + exp(-((r-c*t)/lambda).^2).*Ecr.*(ampC(:,i1)./r.^d);
end

end

