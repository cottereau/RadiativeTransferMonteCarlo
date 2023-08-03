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

% test
if numel(y)>1
    error('the points should be in plane (x,z) and y should be a scalar')
end

% constants
minR = eps;
Nt = length(t);

% initialization
Ec = zeros(numel(x),Nt);
if d==2
    y=0; 
    y0 = zeros(size(x0));
end

% loop on sources
for i1 = 1:lenth(x0)
    r = sqrt((x(:)-x0(1)).^2+(y-y0(i1)).^2+(z(:)-z0(1)).^2);
    r(r<minR)=minR;
    Ec = Ec + exp(-((r-c*t)/lambda).^2).*(Ecr./r.^d);
end
