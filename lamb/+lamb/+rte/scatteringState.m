function S = scatteringState(iface, inc, p, omega, z, opts)
%SCATTERINGSTATE  Field seen by a scatterer near an interface (distorted-wave Born state).
%   S = lamb.rte.scatteringState(iface, inc, p, omega, z)
%   S = lamb.rte.scatteringState(iface, inc, p, omega, z, 'Mode', 'incident')
%   S = lamb.rte.scatteringState(iface, 'IW1', [], omega, z)
%
%   inc   incident body wave of UNIT displacement amplitude:
%           'P1' 'SV1' 'SH1'  from the lower solid (z > 0), travelling up toward z = 0;
%           'Pf'              from the fluid (z < 0), travelling down;
%           'P2' 'SV2' 'SH2'  from an upper solid (z < 0), travelling down;
%         or 'IWj', the j-th real interface-wave mode (lamb.findPoles order), normalised to
%         unit power per unit wavefront width.
%   p     horizontal slowness (ignored for 'IWj'). The frame is the incident's local frame:
%         x along the horizontal propagation direction, y = t, z down (lambConventionMap).
%   z     depths (vector); z > 0 lower medium, z < 0 upper medium.
%   Mode  'state' (default): incident + reflected + transmitted waves (half-space scattering
%         state, the DWBA field); 'incident': the incident plane wave only (the convention of
%         Xu et al. 2021 B.3-B.5 and of scholteModel/stoneleyModel).
%   Away  true: a single plane wave travelling AWAY from the interface in the same medium
%         (the field of a particle emitted inside the slab on its way out); implies 'incident'.
%
%   S.u    3-by-Nz complex displacement at (x = 0, z), horizontal phase e^{i k x} factored out
%   S.eps  3-by-3-by-Nz complex strain
%   S.I    energy flux of the incident wave per unit area normal to it, 0.5 rho c omega^2;
%          for 'IWj', the power per unit wavefront width of the mode (= 1)
%   S.k    horizontal wavenumber
%
%   Polarisations (unit displacement): P along the propagation direction, SV = dir x t,
%   SH = t, as in lambConventionMap. Potentials: A = -1i/k_alpha, B = -1i/k_beta, C = 1.
%
%   References:
%     - Xu, Z., Margerin, L. & Mikesell, T.D. (2021), Monte Carlo simulations of coupled
%       body- and Rayleigh-wave multiple scattering in elastic media, Geophys. J. Int. 228,
%       1213-1236. doi:10.1093/gji/ggab370
arguments
    iface (1,1) lamb.Interface
    inc char
    p double
    omega (1,1) double
    z double
    opts.Mode char {mustBeMember(opts.Mode, {'state', 'incident'})} = 'state'
    opts.Away (1,1) logical = false
end
z = z(:).';
if strncmp(inc, 'IW', 2)
    S = modeState(iface, str2double(inc(3:end)), omega, z);
    return
end
if any(strcmp(inc, {'P2', 'SV2', 'SH2'}))                % upper solid: mirror the problem
    flip = lamb.Interface(iface.upper, iface.lower);
    S = lamb.rte.scatteringState(flip, [inc(1:end-1) '1'], p, omega, -z, 'Mode', opts.Mode, ...
                                 'Away', opts.Away);
    Mz = diag([1 1 -1]);
    S.u = Mz*S.u;
    S.eps = pagemtimes(pagemtimes(Mz, S.eps), Mz);
    return
end

k = omega*p;
L = iface.lower;  U = iface.upper;
S.k = k;  S.u = zeros(3, numel(z));  S.eps = zeros(3, 3, numel(z));
lo = z >= 0;  up = ~lo;
if strcmp(inc, 'Pf')
    kf = U.wavenumber(omega, 'P');
    amp = struct('Pf0', -1i/kf);
    S.I = 0.5*U.rho*U.alpha*omega^2;
    if opts.Away                                          % fluid P travelling up, away
        S = addWave(S, k, U.verticalWavenumber(k, omega, 'P'), +1, 'P', -1i/kf, z, up);
        return
    end
    S = addWave(S, k, U.verticalWavenumber(k, omega, 'P'), -1, 'P', amp.Pf0, z, up);   % incident
    if strcmp(opts.Mode, 'state')
        x = iface.solve(amp, k, omega);
        S = addWave(S, k, U.verticalWavenumber(k, omega, 'P'), +1, 'P', x.Af, z, up);
        S = addWave(S, k, L.verticalWavenumber(k, omega, 'P'), -1, 'P',  x.A1, z, lo);
        S = addWave(S, k, L.verticalWavenumber(k, omega, 'S'), -1, 'SV', x.B1, z, lo);
    end
    return
end

% incident from the lower solid, travelling up (e^{+nu z})
na = L.verticalWavenumber(k, omega, 'P');  nb = L.verticalWavenumber(k, omega, 'S');
amp = struct('P0', 0, 'S0', 0, 'C0', 0);
switch inc
    case 'P1',  amp.P0 = -1i/L.wavenumber(omega, 'P');  S.I = 0.5*L.rho*L.alpha*omega^2;
    case 'SV1', amp.S0 = -1i/L.wavenumber(omega, 'S');  S.I = 0.5*L.rho*L.beta*omega^2;
    case 'SH1', amp.C0 = 1;                             S.I = 0.5*L.rho*L.beta*omega^2;
    otherwise,  error('lamb:scatteringState:inc', 'Unknown incident wave %s.', inc);
end
if opts.Away                                              % travelling down, away (e^{-nu z})
    S = addWave(S, k, na, -1, 'P',  amp.P0, z, lo);
    S = addWave(S, k, nb, -1, 'SV', amp.S0, z, lo);
    S = addWave(S, k, nb, -1, 'SH', amp.C0, z, lo);
    return
end
S = addWave(S, k, na, +1, 'P',  amp.P0, z, lo);
S = addWave(S, k, nb, +1, 'SV', amp.S0, z, lo);
S = addWave(S, k, nb, +1, 'SH', amp.C0, z, lo);
if strcmp(opts.Mode, 'incident'), return, end
x = iface.solve(amp, k, omega);
S = addWave(S, k, na, -1, 'P',  x.A1, z, lo);            % reflected
S = addWave(S, k, nb, -1, 'SV', x.B1, z, lo);
S = addWave(S, k, nb, -1, 'SH', x.C1, z, lo);
switch iface.kind
    case 'solid-solid'                                   % transmitted, e^{+nu z} for z < 0
        S = addWave(S, k, U.verticalWavenumber(k, omega, 'P'), +1, 'P',  x.A2, z, up);
        S = addWave(S, k, U.verticalWavenumber(k, omega, 'S'), +1, 'SV', x.B2, z, up);
        S = addWave(S, k, U.verticalWavenumber(k, omega, 'S'), +1, 'SH', x.C2, z, up);
    case 'solid-fluid'
        S = addWave(S, k, U.verticalWavenumber(k, omega, 'P'), +1, 'P',  x.Af, z, up);
end
end

% ------------------------------------------------------------------------------------------
function S = modeState(iface, j, omega, z)
% Interface-wave eigenstate from the null vector of M(k_p), unit power per unit width.
poles = lamb.findPoles(iface, omega);
kp = poles(j).k;
r = lamb.poleNullVectors(iface, kp, omega);
names = iface.unknowns();  res = struct();
for i = 1:numel(names), res.(names{i}) = r(i); end
r = r/sqrt(lamb.energy.interfaceFlux2D(iface, kp, omega, res));
L = iface.lower;  U = iface.upper;
S.k = kp;  S.I = 1;  S.u = zeros(3, numel(z));  S.eps = zeros(3, 3, numel(z));
lo = z >= 0;  up = ~lo;
S = addWave(S, kp, L.verticalWavenumber(kp, omega, 'P'), -1, 'P',  r(1), z, lo);
S = addWave(S, kp, L.verticalWavenumber(kp, omega, 'S'), -1, 'SV', r(2), z, lo);
switch iface.kind
    case 'solid-solid'
        S = addWave(S, kp, U.verticalWavenumber(kp, omega, 'P'), +1, 'P',  r(3), z, up);
        S = addWave(S, kp, U.verticalWavenumber(kp, omega, 'S'), +1, 'SV', r(4), z, up);
    case 'solid-fluid'
        S = addWave(S, kp, U.verticalWavenumber(kp, omega, 'P'), +1, 'P',  r(3), z, up);
end
end

function S = addWave(S, k, nu, s, type, a, z, mask)
% Add the plane wave (potential or SH amplitude a) proportional to e^{i k x + s nu z}, at
% depths z(mask): displacement u and gradient G_ij = du_i/dx_j = u_i kappa_j.
if ~any(mask) || a == 0, return, end
e = exp(s*nu*z(mask));
kap = [1i*k; 0; s*nu];
switch type
    case 'P',  pol = [1i*k; 0; s*nu];                    % grad(phi)
    case 'SV', pol = [-s*nu; 0; 1i*k];                   % curl(psi*y)
    case 'SH', pol = [0; 1; 0];
end
u = a*pol*e;                                             % 3-by-n
S.u(:, mask) = S.u(:, mask) + u;
G = reshape(u, 3, 1, []).*reshape(kap, 1, 3);            % 3-by-3-by-n
S.eps(:, :, mask) = S.eps(:, :, mask) + (G + pagetranspose(G))/2;
end
