function varargout = lambConventionMap(varargin)
%LAMBCONVENTIONMAP  MC particle state <-> local frame of the +lamb f-k formulation.
%
%   L = lambConventionMap(dir, perp, nrm)       MC state -> Lamb frame
%   [dir, perp] = lambConventionMap(L)          Lamb frame -> MC state
%
%   dir, perp : N-by-3 unit propagation and polarisation vectors (perp . dir = 0), as in
%               P.dir / P.perp of the MC engines. P particles carry a perp too; it is mapped
%               the same way, so the round trip is exact for every particle.
%   nrm       : 1-by-3 or N-by-3 unit interface normal pointing INTO Lamb medium 1, the
%               half-space that is "z > 0" in +lamb (the solid of a solid-fluid pair; for a
%               solid-solid pair, the side of the incident particle).
%
%   Lamb frame (LambFK_Formulation.tex, Secs. 1 and 7):
%     zhat = nrm,   lhat = horizontal part of dir (normalised),   that = zhat x lhat,
%     dir  = sinT*lhat + cosT*zhat,  cosT = dir . zhat  (cosT < 0: toward medium 2, the
%     "up-going" P0/S0/C0 of +lamb; cosT > 0: away from the interface into medium 1).
%   (lhat, that, zhat) is right-handed like (x, y, z), so the 2D out-of-plane axis y of
%   u = grad(phi) + curl(psi*y) is that.
%
%   Unit-displacement polarisations of a plane wave e^{i k d.x}:
%     P  : u = dir                    potential  A = aP /(1i*k_alpha)
%     SV : u = eSV = dir x that       potential  B = aSV/(1i*k_beta)   (u = curl(psi*that))
%     SH : u = that                   displacement C = aSH
%   with aSV = perp . eSV and aSH = perp . that, so aSV^2 + aSH^2 = 1.
%
%   Output struct L (all N-by-3 or N-by-1): .zhat .lhat .that .eSV .sinT .cosT .aSV .aSH
%
%   Normal incidence (horizontal part of dir exactly zero): lhat is taken along perp, so all
%   S energy is SV (aSH = 0). SV and SH are degenerate there; any tangent is correct.
%
%   See also lamb.incidentAmplitudes.

if nargin == 1                                   % Lamb frame -> MC state
    L = varargin{1};
    dir  = L.sinT.*L.lhat + L.cosT.*L.zhat;
    perp = L.aSV.*L.eSV + L.aSH.*L.that;
    varargout = {dir, perp};
    return
end

[dir, perp, nrm] = varargin{1:3};
N = size(dir, 1);
z = nrm.*ones(N, 1);                             % broadcast a single normal
cosT = sum(dir.*z, 2);
h    = dir - cosT.*z;                            % horizontal part of dir
sinT = vecnorm(h, 2, 2);
l    = h./sinT;
at0  = sinT == 0;                                % exact normal incidence
if any(at0)
    pt = perp(at0, :) - sum(perp(at0, :).*z(at0, :), 2).*z(at0, :);
    l(at0, :) = pt./vecnorm(pt, 2, 2);
end
t   = cross(z, l, 2);
eSV = cross(dir, t, 2);

L = struct('zhat', z, 'lhat', l, 'that', t, 'eSV', eSV, 'sinT', sinT, 'cosT', cosT, ...
           'aSV', sum(perp.*eSV, 2), 'aSH', sum(perp.*t, 2));
varargout = {L};
end
