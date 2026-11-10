function out = sourceVector(varargin)
%SOURCEVECTOR  Pack / unpack the 9-component source vector.
%   s = lamb.rte.sourceVector(src)          % from a lamb.PointSource
%   s = lamb.rte.sourceVector(F, M)         % F 3-by-1, M 3-by-3 symmetric (either may be [])
%   B = lamb.rte.sourceVector(s)            % unpack: B.F (3-by-1), B.M (3-by-3 symmetric)
%
%   s = [F_x F_y F_z  M_xx M_yy M_zz M_xy M_xz M_yz].'  The off-diagonal basis elements are
%   E_pq + E_qp, so the coefficient of E_xy + E_yx is the matrix entry M_xy. Every radiated
%   amplitude is linear in s; energies are Hermitian forms s' * Q * s.
if nargin == 1 && isa(varargin{1}, 'lamb.PointSource')
    src = varargin{1};
    switch src.kind
        case 'force',  out = pack(src.F, []);
        case 'moment', out = pack([], src.M);
        otherwise,     out = pack(src.F, src.M);
    end
elseif nargin == 2
    out = pack(varargin{1}, varargin{2});
else
    s = varargin{1};
    out.F = s(1:3);  out.F = out.F(:);
    out.M = [s(4) s(7) s(8); s(7) s(5) s(9); s(8) s(9) s(6)];
end
end

function s = pack(F, M)
if isempty(F), F = zeros(3, 1); end
if isempty(M), M = zeros(3); end
if norm(M - M.', 'fro') > 1e-12*max(1, norm(M, 'fro'))
    error('lamb:sourceVector:asym', 'M must be symmetric (the torque part is not in the basis).');
end
s = [F(:); M(1,1); M(2,2); M(3,3); M(1,2); M(1,3); M(2,3)];
end
