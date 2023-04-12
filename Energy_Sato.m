function E = Energy_Sato(d,Sigma,v,r,t)

% d     : dimension of the problem
% Sigma : total scattering cross-section
% v     : velocity (acoustic, shear elastic)
% r     : source-receiver distance
% t     : time vector

E = zeros(size(t));

if d==2

    E = Sigma*exp(Sigma*sqrt(v^2*t.^2-r^2)-Sigma*v*t).*heaviside(t-r/v) ...
        ./ (2*pi*sqrt(v^2*t.^2-r^2));
    
    ind = (t==r/v);
    if ~isempty(E(ind))
        E(ind) = E(ind) + exp(-Sigma*v*t(ind))/(2*pi*v*r);
    end

elseif d==3
    pause 2 % to be changed
else
    disp('error : input dimension should be either 2 or 3')
end
