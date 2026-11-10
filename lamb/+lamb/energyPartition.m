function EP = energyPartition(iface, src, omega, opts)
%ENERGYPARTITION  Radiated-power partition between body waves and interface waves.
%   EP = lamb.energyPartition(iface, src, omega)
%   EP = lamb.energyPartition(iface, src, omega, 'NGauss', 64, 'NPhi', 8, ...)
%
%   Evaluated at REAL omega (lossless). Real poles (Rayleigh / Scholte / Stoneley) are
%   removed analytically by residues; the body-wave windows contain no real poles. (The
%   complex frequency omega + 1i*eps is for wavefield synthesis only: it dissipates energy
%   and moves the poles off the axis, which blurs the interface/body split.)
%
%   Options
%     NGauss        Gauss nodes per quadrature segment (default 64)
%     NPhi          azimuth nodes, >= 5 is exact for forces and moment tensors (default 8)
%     CheckBalance  compute the source work and the relative error (default true)
%
%   Output EP
%     .omega
%     .W.body       struct [W]: .P1 .SV1 .SH1 (+ .P2 .SV2 .SH2 | .Pf)
%     .W.interface  struct array per real pole: .type .k .c .W (flux) .Wsrc (source side)
%     .W.leaky      struct array per near-real leaky pole: .type .k .leakLength .W
%                   (W is the part of W.body inside the Lorentzian peak: counted in the body
%                   total, NOT added again; it stays trapped over leakLength)
%     .W.bodyTotal, .W.ifTotal, .W.total (= bodyTotal + ifTotal)
%     .pct          same layout in percent of W.total: .body.(name), .interface (vector),
%                   .leaky (vector), .bodyTotal, .ifTotal, .trapped (= ifTotal + sum(leaky))
%     .check        .Wsrc (struct from lamb.energy.sourcePower), .relErr,
%                   .ifRelErr (flux vs source-side value, per real pole)
arguments
    iface (1,1) lamb.Interface
    src   (1,1) lamb.PointSource
    omega (1,1) double {mustBeReal, mustBePositive}
    opts.NGauss (1,1) double = 64
    opts.NPhi   (1,1) double = 8
    opts.CheckBalance (1,1) logical = true
end
[poles, leaky] = lamb.findPoles(iface, omega);
[Wb, W0, Wleak] = lamb.energy.bodyWaves(iface, src, omega, poles, leaky, opts);
Wi = lamb.energy.interfaceWaves(iface, src, poles, omega, opts);

EP.omega = omega;
EP.W.body      = Wb;
EP.W.interface = Wi;
EP.W.leaky     = struct('type', {}, 'k', {}, 'leakLength', {}, 'W', {});
for j = 1:numel(leaky)
    EP.W.leaky(j) = struct('type', leaky(j).type, 'k', leaky(j).k, ...
                           'leakLength', leaky(j).leakLength, 'W', Wleak(j));
end
EP.W.bodyTotal = sum(cell2mat(struct2cell(Wb)));
EP.W.ifTotal   = sum([Wi.W]);
EP.W.total     = EP.W.bodyTotal + EP.W.ifTotal;

pc = @(v) 100*v/EP.W.total;
EP.pct.body      = structfun(pc, Wb, 'UniformOutput', false);
EP.pct.interface = pc([Wi.W]);
EP.pct.leaky     = pc(Wleak);
EP.pct.bodyTotal = pc(EP.W.bodyTotal);
EP.pct.ifTotal   = pc(EP.W.ifTotal);
EP.pct.trapped   = pc(EP.W.ifTotal + sum(Wleak));

if opts.CheckBalance
    Ws = lamb.energy.sourcePower(iface, src, omega, W0, Wi, poles, leaky, opts);
    EP.check.Wsrc     = Ws;
    EP.check.relErr   = abs(Ws.total - EP.W.total)/EP.W.total;
    EP.check.ifRelErr = abs([Wi.W] - [Wi.Wsrc])./abs([Wi.W]);
end
end
