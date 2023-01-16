function P = position2angle(d,P,lradius)
% P = POSITION2ANGLE(d,P,lradius) computes the angle P.theta and radius P.r
% corresponding to the cylindrical coordinates of point (X,Y) in d=2, or 
% the angles P.theta and P.phi and radius P.r corresponding to the 
% spherical coordinates of point (X,Y,Z) in d=3.
% 
% When lradius=true, functions CART2POL is used in 2D, and CART2SPH in 3D.
%
% When lradius=false, the routine expects to find the radius already
% computed in P.r.

% 2D
if d==2
    if lradius
        P.theta = acos(P.x./P.r);
        as = asin(P.y./P.r);
        P.theta(as<0) = -P.theta(as<0);
    else
        [P.theta,P.r] = cart2pol(P.x,P.y);
    end

% 3D to be done: computing the angle in 3D
elseif d==3
    warning('POSITION2ANGLE: routine not implemented in 3D yet')
else
    warning('POSITION2ANGLE: routine only works for 2D and 3D')
end