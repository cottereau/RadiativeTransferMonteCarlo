function [t, u, info] = synthetics(iface, src, rcv, opts)
%SYNTHETICS  Time-domain displacement seismograms with the complex-frequency method.
%   [t, u] = lamb.wavefield.synthetics(iface, src, rcv, 'T', T, 'dt', dt)
%   [t, u, info] = lamb.wavefield.synthetics(..., 'Eps', eps, 'L', L, 'kMax', kMax)
%
%   T, dt  time window and sampling [s]; Nt = 2^nextpow2(T/dt)
%   Eps    imaginary frequency (default pi/T): omega_c = omega + 1i*Eps. The factor
%          exp(Eps*t) applied after synthesis undoes the damping exactly (to rounding), so the
%          seismograms are those of the undamped problem.
%   L      DWN periodicity (default 1.1*(r_max + c_max*T))
%   kMax, Direct   passed to lamb.wavefield.dwn
%   Frequencies where |wavelet| < 1e-10*max|wavelet| are skipped (spectrum set to 0).
%   Parts, SplitSlowness   field split, passed to lamb.wavefield.dwn
%   SourceFun  frequency-dependent source vector, passed to lamb.wavefield.dwn (src gives zs and
%          the wavelet, which still multiplies the spectra)
%   u      Nt-by-3-by-Nr (u_x, u_y, u_z), convolved with src.wavelet (Nt-by-3-by-Nr-by-numel(Parts)
%          with Parts)
%   info   .omegaC .L .eps
arguments
    iface (1,1) lamb.Interface
    src   (1,1) lamb.PointSource
    rcv   struct
    opts.T    (1,1) double {mustBePositive}
    opts.dt   (1,1) double {mustBePositive}
    opts.Eps  (1,1) double = NaN
    opts.L    (1,1) double = NaN
    opts.kMax (1,1) double = NaN
    opts.Direct (1,1) logical = true
    opts.Parts (1,:) string = strings(1, 0)
    opts.SplitSlowness (1,1) double = NaN
    opts.SourceFun = []
end
Nt = 2^nextpow2(ceil(opts.T/opts.dt));
T  = Nt*opts.dt;
if isnan(opts.Eps), opts.Eps = pi/T; end
if isnan(opts.L)
    c = [iface.lower.alpha, iface.lower.beta];
    if ~isempty(iface.upper), c = [c, iface.upper.alpha]; end
    opts.L = 1.1*(max(hypot(rcv.x(:), rcv.y(:))) + max(c)*T);
end
omegaC = (0:Nt/2).'*(2*pi/T) + 1i*opts.Eps;
S    = src.wavelet(omegaC);
keep = abs(S) > 1e-10*max(abs(S));                          % skip frequencies outside the band
U = zeros(numel(rcv.x), 3, numel(omegaC), max(1, numel(opts.Parts)));
U(:, :, keep, :) = lamb.wavefield.dwn(iface, src, rcv, omegaC(keep), 'L', opts.L, ...
                                   'kMax', opts.kMax, 'Direct', opts.Direct, ...
                                   'Parts', opts.Parts, 'SplitSlowness', opts.SplitSlowness, ...
                                   'SourceFun', opts.SourceFun);
U = U.*reshape(S, 1, 1, []);
[t, u] = lamb.wavefield.freq2time(permute(U, [3 2 1 4]), opts.dt, opts.Eps);   % Nt-by-3-by-Nr
info = struct('omegaC', omegaC, 'L', opts.L, 'eps', opts.Eps);
end
