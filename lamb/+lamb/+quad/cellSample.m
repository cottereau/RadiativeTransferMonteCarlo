function k = cellSample(cells, u)
%CELLSAMPLE  Map uniform numbers u into k inside quadrature cells (lamb.quad.windowNodes).
%   k = lamb.quad.cellSample(cells, u)   cells: M-by-5 rows [tLo tHi ka kb type], u: M-by-1
%   in [0, 1). t = tLo + u (tHi - tLo) is uniform in the segment's mapped variable, where the
%   integrand is smooth (the map removes the branch-point singularities), then mapped to k.
t  = cells(:, 1) + u.*(cells(:, 2) - cells(:, 1));
ka = cells(:, 3);  kb = cells(:, 4);
k  = sqrt(ka.^2 + (kb.^2 - ka.^2).*sin(t).^2);                 % type 1
lin = cells(:, 5) == 2;
k(lin) = (ka(lin) + kb(lin))/2 + (kb(lin) - ka(lin))/2.*t(lin);
end
