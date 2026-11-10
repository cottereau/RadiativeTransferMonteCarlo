function B = basisSources(zs)
%BASISSOURCES  The 9 basis point sources of the source vector s (lamb.rte.sourceVector).
%   B = lamb.rte.basisSources(zs)  -> 1-by-9 cell of lamb.PointSource at depth zs:
%   forces e_x, e_y, e_z, then moments E_xx, E_yy, E_zz, E_xy+E_yx, E_xz+E_zx, E_yz+E_zy.
B = cell(1, 9);
for i = 1:9
    e = zeros(9, 1);  e(i) = 1;
    u = lamb.rte.sourceVector(e);
    if i <= 3
        B{i} = lamb.PointSource('force', u.F, zs);
    else
        B{i} = lamb.PointSource('moment', u.M, zs);
    end
end
end
