function [poles, leaky] = findPoles(iface, omega, opts)
%FINDPOLES  Interface-wave poles at real omega.
%   [poles, leaky] = lamb.findPoles(iface, omega)
%
%   poles : real poles k_p > max(branch points) on the physical sheet (trapped waves).
%           Found by scanning the real dispersion function for sign changes, then fzero.
%           Fields .k .c .type ('Rayleigh' | 'Scholte' | 'Stoneley')
%   leaky : near-real LEAKY poles, k_L = a + 1i*b with 0 < b/a < opts.LeakyTol, lying inside
%           a propagating window (between the solid's k_beta and the largest branch point),
%           e.g. the leaky Rayleigh pole of a solid under a slow fluid. On the real axis they
%           appear as Lorentzian peaks of half-width b, which the quadrature must resolve
%           (see lamb.quad.windowNodes). Their energy is body-wave energy at infinite range,
%           but it stays trapped over the leakage length 1/(2b).
%           Fields .k (complex) .a .b .leakLength .type
%
%   Leaky poles are located with real-axis Newton steps k <- a - D(a)/D'(a), a <- Re(k):
%   the real-axis dispersion function is analytic inside each window, so the linear model
%   gives Im(k_L) without evaluating any function at complex k.
arguments
    iface (1,1) lamb.Interface
    omega (1,1) double {mustBeReal, mustBePositive}
    opts.kMaxFactor (1,1) double = 5
    opts.NScan      (1,1) double = 4000
    opts.NSeed      (1,1) double = 24
    opts.LeakyTol   (1,1) double = 1e-3
end
bp   = iface.branchPoints(omega);
name = iface.interfaceWave();
leakName = struct('free','leaky Rayleigh','solid_fluid','leaky Rayleigh', ...
                  'solid_solid','leaky Stoneley').(strrep(iface.kind, '-', '_'));
D    = @(k) iface.dispersion(k, omega);

% --- real (trapped) poles -----------------------------------------------------------
k0 = bp(end);
ks = k0*(1 + logspace(-10, log10(opts.kMaxFactor - 1), opts.NScan)).';
Ds = real(D(ks));
idx = find(sign(Ds(1:end-1)).*sign(Ds(2:end)) < 0);
poles = struct('k', {}, 'c', {}, 'type', {});
for j = idx.'
    kp = fzero(@(k) real(D(k)), ks([j j+1]));
    poles(end+1) = struct('k', kp, 'c', omega/kp, 'type', name); %#ok<AGROW>
end

% --- near-real leaky poles inside windows beyond the lower solid's k_beta ---------------
leaky = struct('k', {}, 'a', {}, 'b', {}, 'leakLength', {}, 'type', {});
kb1 = iface.lower.wavenumber(omega, 'S');
segs = bp(bp >= kb1);
for s = 1:numel(segs) - 1
    lo = segs(s);  hi = segs(s+1);
    for a = lo*(hi/lo).^(((1:opts.NSeed) - 0.5)/opts.NSeed)      % geometric seeds
        ok = false;
        for it = 1:60
            h  = 1e-7*a;
            kL = a - D(a)/((D(a + h) - D(a - h))/(2*h));
            if real(kL) <= lo || real(kL) >= hi, break, end
            if abs(real(kL) - a) < 1e-15*a, ok = true; break, end
            a = real(kL);
        end
        b = abs(imag(kL));
        if ok && b < opts.LeakyTol*a && b > 0 && ...
                (isempty(leaky) || all(abs([leaky.a] - a) > 1e-9*a))
            leaky(end+1) = struct('k', kL, 'a', a, 'b', b, 'leakLength', 1/(2*b), ...
                                  'type', leakName); %#ok<AGROW>
        end
    end
end
end
