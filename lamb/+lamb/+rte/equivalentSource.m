function s = equivalentSource(S, het, omega)
%EQUIVALENTSOURCE  Point source equivalent to a small heterogeneity in a given field (Born).
%   s = lamb.rte.equivalentSource(S, het, omega)
%
%   S    field at the scatterer (lamb.rte.scatteringState): S.u (3-by-Nz), S.eps (3-by-3-by-Nz)
%   het  struct with .V (volume) and the ABSOLUTE contrasts .drho, .dlam, .dmu (scalars or
%        1-by-Nz); missing fields are zero
%   s    9-by-Nz source vectors (lamb.rte.sourceVector order):
%            F = omega^2 drho V u,     M = -V (dlam tr(eps) I + 2 dmu eps),
%        in the convention u_i = M_pq dG_ip/dxi_q of lamb.PointSource, i.e. the Born
%        scattered field int [omega^2 drho u_i G_ni - dC_ijkl eps_kl dG_ni/dx'_j] dV.
V = het.V;
g = @(f) getf(het, f);
drho = g('drho');  dlam = g('dlam');  dmu = g('dmu');
E  = S.eps;
tr = squeeze(E(1,1,:) + E(2,2,:) + E(3,3,:)).';
F  = omega^2*V*drho.*S.u;
M  = @(i, j) -V*(dlam.*tr*(i == j) + 2*dmu.*squeeze(E(i,j,:)).');
s  = [F; M(1,1); M(2,2); M(3,3); M(1,2); M(1,3); M(2,3)];
end

function v = getf(het, f)
if isfield(het, f), v = het.(f); else, v = 0; end
end
