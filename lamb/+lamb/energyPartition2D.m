function EP = energyPartition2D(iface, src, omega, opts)
%ENERGYPARTITION2D  Radiated-power partition of a LINE source (plane strain, 2D P-SV).
%   EP = lamb.energyPartition2D(iface, src, omega)
%   EP = lamb.energyPartition2D(iface, src, omega, 'NGauss', 64)
%
%   The source src (lamb.PointSource) is extended uniformly along y: a line source at depth
%   src.zs, the 2D problem of LambFK_Formulation.tex. Its 3D transform is 2 pi delta(k_y)
%   times the point source's, so every amplitude is the 3D one at k_y = 0, i.e. at the two
%   azimuths phi = 0 (k_x = +k) and phi = pi (k_x = -k), and the powers are per unit length
%   along y:
%     body       W_m = (1/2pi) int_0^{k_m} E_m(k) (|D_m(k,0)|^2 + |D_m(k,pi)|^2) dk
%                (3D: (1/4pi^2) int int E_m |D_m|^2 k dk dphi)
%     interface  W_j = W2D(Res(phi = 0)) + W2D(Res(phi = pi))  (lamb.energy.interfaceFlux2D;
%                the pole term of (1/2pi) int dk_x is i Res e^{i k_p |x|}, one per side)
%   Real poles only; near-real leaky poles are counted in the body powers (no Wleak split).
%   For an in-plane source (F_y = 0, M_xy = M_yz = 0) the SH powers are zero.
%
%   Output EP: .omega, .W.body (struct, as lamb.energyPartition), .W.interface (.type .k
%   .c .W .Wphi = [toward +x, toward -x]), .nodes.(name) (.k quadrature nodes, .dW N-by-2
%   power per node for phi = 0 (+x) and pi (-x): the emission pattern for a 2D Monte Carlo),
%   .W.bodyTotal, .W.ifTotal, .W.total, .pct.body (struct), .pct.interface, .pct.bodyTotal,
%   .pct.ifTotal.
%
%   Reference (Poisson solid, vertical line force at h = 0): P 21.94 %, S 29.13 %,
%   Rayleigh 48.93 % (independent Python evaluation).
%
%   See also LAMB.ENERGYPARTITION, LAMB.ENERGY.INTERFACEFLUX2D.
arguments
    iface (1,1) lamb.Interface
    src   (1,1) lamb.PointSource
    omega (1,1) double {mustBeReal, mustBePositive}
    opts.NGauss (1,1) double = 64
end
phi = [0 pi];
[poles, leaky] = lamb.findPoles(iface, omega);
bp    = iface.branchPoints(omega);
peaks = lamb.quad.peaksFor(iface, omega, poles, leaky);
for w = lamb.energy.outgoingWaves(iface, src.zs >= 0)
    km = w.med.wavenumber(omega, w.type);
    [k, wk] = lamb.quad.windowNodes([0, bp(bp < km), km], opts.NGauss, peaks);
    amp = lamb.incidentAmplitudes(src, iface, k, phi, omega, true);
    x   = iface.solve(amp, k, omega);
    D   = x.(w.field);
    if ~isempty(w.away), D = D + amp.(w.away); end
    dW = wk.*lamb.energy.fluxDensity(w.med, k, omega, w.type).*abs(D).^2/(2*pi);   % N-by-2 (phi = 0, pi)
    Wb.(w.name) = sum(dW, 'all');
    nodes.(w.name) = struct('k', k, 'dW', dW);
end
Wi = struct('type', {}, 'k', {}, 'c', {}, 'W', {}, 'Wphi', {});
for p = poles
    res = lamb.poleResidues(iface, src, p, phi, omega);
    Wphi = lamb.energy.interfaceFlux2D(iface, p.k, omega, res);              % toward +x, -x
    Wi(end+1) = struct('type', p.type, 'k', p.k, 'c', p.c, 'W', sum(Wphi), 'Wphi', Wphi);   %#ok<AGROW>
end

EP.omega = omega;
EP.W.body      = Wb;
EP.W.interface = Wi;
EP.nodes       = nodes;
EP.W.bodyTotal = sum(cell2mat(struct2cell(Wb)));
EP.W.ifTotal   = sum([Wi.W]);
EP.W.total     = EP.W.bodyTotal + EP.W.ifTotal;
pc = @(v) 100*v/EP.W.total;
EP.pct.body      = structfun(pc, Wb, 'UniformOutput', false);
EP.pct.interface = pc([Wi.W]);
EP.pct.bodyTotal = pc(EP.W.bodyTotal);
EP.pct.ifTotal   = pc(EP.W.ifTotal);
end
