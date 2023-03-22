function [corners,x,z] = computBox(geom,pos,dx)
if isfield(geom,'planeX')
    x0 = geom.planeX(:,1);
else
    x0 = inf;
end
x = boxBounds(x0,geom.movieBox(1),pos(1));
if isfield(geom,'planeZ')
    z0 = geom.planeZ(:,1);
else
    z0 = inf;
end
z = boxBounds(z0,geom.movieBox(2),pos(3));
corners = [x(1) 0 z(1); 
           x(2) 0 z(1); 
           x(2) 0 z(2); 
           x(1) 0 z(2)];

if nargin>2
    x = linspace(x(1),x(2),ceil((x(2)-x(1))/dx));
    z = linspace(z(1),z(2),ceil((z(2)-z(1))/dx));
end
end

%%%%%%%%%%%%%%%%%%%
function x = boxBounds(x0,L,pos)
if length(x0)==1
    if isinf(x0)
        x = [-L/2 L/2];
    elseif x0>pos
        x = [x0-L x0];
    else
        x = [x0 x0+L];
    end
elseif length(x0)==2
    x = sort(x0)';
end
end
