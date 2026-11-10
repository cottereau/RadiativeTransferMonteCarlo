function U = dwn(iface, src, rcv, omegaC, opts)
%DWN  Displacement spectra at receivers by the discrete wavenumber method (Bouchon).
%   U = lamb.wavefield.dwn(iface, src, rcv, omegaC, 'L', L)
%   U = lamb.wavefield.dwn(..., 'kMax', kMax, 'NPhi', 8)
%
%   rcv     struct with column vectors .x .y .z [m] (z positive down; z >= 0 is the solid)
%   omegaC  complex frequencies omega + 1i*eps (eps > 0), vector. The imaginary part moves
%           the real poles off the k axis, so the discrete k sum is regular.
%   L       source periodicity [m]: k_j = j*2*pi/L, j = 1..Nk. Periodic images arrive after
%           (L - r)/c_max; choose L > r_max + c_max*T (lamb.wavefield.synthetics does this).
%   kMax    largest wavenumber (default: 1.5*max|omega|/c_min + 25/d_min, where d_min is the
%           smallest vertical distance governing the evanescent decay; set it explicitly
%           when a receiver is at the source depth).
%   Direct  include the whole-space direct field (default true); false = scattered field only
%           (reflected/transmitted waves), e.g. for receivers at the source depth.
%   Parts   optional split of the field, e.g. ["P" "SV" "SH" "IW"]; U gets a 4th dimension
%           (one page per part, the pages sum to the total). P, SV, SH: the P, SV, SH
%           potentials (direct + reflected/transmitted; fluid: P) at k <= real(omega)*SplitSlowness,
%           SH at all k; IW: the P and SV potentials at k > real(omega)*SplitSlowness, i.e.
%           evanescent in every medium, where the interface-wave pole lies (plus the near field
%           of the source). Give SplitSlowness between 1/c_body,min and 1/c_IW.
%   SourceFun  optional @(omega) -> 9-vector (lamb.rte.sourceVector order, complex): a source whose
%           force/moment depend on frequency (e.g. the Born source of a scatterer in an incident
%           wave); at each omega the source is lamb.PointSource('general', SourceFun(omega), src.zs)
%   U       Nr-by-3-by-Nw complex spectra (u_x, u_y, u_z), NOT multiplied by the wavelet
%           (Nr-by-3-by-Nw-by-numel(Parts) with Parts).
%
%   Method, per frequency (vectorised over k, phi and receivers):
%     1. amplitudes on the (k, phi) grid: incidentAmplitudes -> Interface.solve (pages)
%     2. fields at each receiver depth in the local frame (u_L, u_T, u_z) -> Cartesian
%     3. FFT over phi -> harmonics n = -3..3 (NPhi >= 8 avoids aliasing)
%     4. u(r,theta,z) = (dk/2pi) sum_n i^n e^{i n theta} sum_j u_n(k_j,z) J_n(k_j r) k_j
arguments
    iface  (1,1) lamb.Interface
    src    (1,1) lamb.PointSource
    rcv    struct
    omegaC double
    opts.L    (1,1) double {mustBePositive}
    opts.kMax (1,1) double = NaN
    opts.NPhi (1,1) double = 8
    opts.Direct (1,1) logical = true
    opts.Parts (1,:) string = strings(1, 0)
    opts.SplitSlowness (1,1) double = NaN
    opts.SourceFun = []
end
np = max(1, numel(opts.Parts));
if np > 1 || ~isempty(opts.Parts)
    assert(all(ismember(opts.Parts, ["P" "SV" "SH" "IW"])), 'lamb:dwn:parts', 'Parts: P, SV, SH, IW.');
    assert(~isnan(opts.SplitSlowness), 'lamb:dwn:split', 'Parts needs SplitSlowness.');
end
x = rcv.x(:);  y = rcv.y(:);  z = rcv.z(:);  Nr = numel(x);
if strcmp(iface.kind, 'free') && any(z < 0)
    error('lamb:dwn:receiver', 'Receivers must be at z >= 0 below a free surface.');
end
r = hypot(x, y);  th = atan2(y, x);
if isnan(opts.kMax)
    c = [iface.lower.alpha, iface.lower.beta];
    if ~isempty(iface.upper), c = [c, iface.upper.alpha, iface.upper.beta]; end
    cmin = 0.85*min(c(c > 0));                               % interface waves are slower
    d = min(abs(z) + abs(src.zs));
    if opts.Direct, d = min([d; abs(z - src.zs)]); end
    if d == 0
        error('lamb:dwn:kMax', 'A receiver is at the source depth or both are at z = 0: give kMax.');
    end
    opts.kMax = 1.5*max(real(omegaC))/cmin + 25/d;
end
dk = 2*pi/opts.L;
k  = (1:ceil(opts.kMax/dk)).'*dk;
[phi, ~] = lamb.quad.azimuthNodes(opts.NPhi);
n  = [0:opts.NPhi/2-1, -opts.NPhi/2:-1];                     % FFT index -> harmonic order
use = abs(n) <= 3;

% Bessel kernels depend on k and r only: build once for all frequencies.
[zu, ~, iz] = unique(z);
J = cell(numel(zu), 1);
for d = 1:numel(zu)
    rr = r(iz == d).';
    for m = find(use)
        J{d}{m} = besselj(n(m), k.*rr);                      % Nk-by-Nr (J_{-n} = (-1)^n J_n)
    end
end

U = zeros(Nr, 3, numel(omegaC), np);
for iw = 1:numel(omegaC)
    w   = omegaC(iw);
    sw  = src;
    if ~isempty(opts.SourceFun), sw = lamb.PointSource('general', opts.SourceFun(w), src.zs); end
    amp = lamb.incidentAmplitudes(sw, iface, k, phi, w);
    xs  = iface.solve(amp, k, w);
    for d = 1:numel(zu)
        idx = find(iz == d);
        for p = 1:np
            if isempty(opts.Parts)
                wt = [1 1 1];  km = true(size(k));
            else
                body = k <= real(w)*opts.SplitSlowness;
                switch opts.Parts(p)
                    case "P",  wt = [1 0 0];  km = body;
                    case "SV", wt = [0 1 0];  km = body;
                    case "SH", wt = [0 0 1];  km = true(size(k));
                    case "IW", wt = [1 1 0];  km = ~body;
                end
            end
            [uL, uT, uzz] = fieldAtDepth(iface, sw, k, phi, w, xs, zu(d), opts.Direct, wt);
            comp = {(uL.*cos(phi) - uT.*sin(phi)).*km, (uL.*sin(phi) + uT.*cos(phi)).*km, uzz.*km};
            for c3 = 1:3
                Un  = fft(comp{c3}, [], 2)/opts.NPhi;        % harmonics (Nk-by-NPhi)
                % Euler-Maclaurin end correction at k = 0: only n = 0 has g'(0) = u_0(0) ~= 0,
                % with u_0(0) from Richardson (u_0 is even in k). Error O(dk^2) -> O(dk^4).
                acc = dk/12*(4*Un(1, 1) - Un(2, 1))/3;
                for m = find(use)
                    acc = acc + (1i^n(m))*exp(1i*n(m)*th(idx).') .* ((k.*Un(:, m)).' * J{d}{m});
                end
                U(idx, c3, iw, p) = dk/(2*pi)*acc.';
            end
        end
    end
end
end

function [uL, uT, uz] = fieldAtDepth(iface, src, k, phi, w, x, z, direct, wt)
% Field in the local frame (l, t, z) at depth z, N-by-Nphi; wt = weights of the P, SV, SH
% potentials ([1 1 1]: total field).
zer = zeros(numel(k), numel(phi));
if z >= 0
    m  = iface.lower;
    na = m.verticalWavenumber(k, w, 'P');  nb = m.verticalWavenumber(k, w, 'S');
    [uL, uz] = psv(k, na, nb, wt(1)*x.A1.*exp(-na*z), wt(2)*x.B1.*exp(-nb*z), -1);
    uT = wt(3)*x.C1.*exp(-nb*z);
    if direct && src.zs >= 0, [uL, uT, uz] = addDirect(uL, uT, uz, m, src, k, phi, w, z, wt); end
else
    m = iface.upper;
    if m.isFluid
        nf = m.verticalWavenumber(k, w, 'P');
        uL = wt(1)*1i*k.*x.Af.*exp(nf*z);  uz = wt(1)*nf.*x.Af.*exp(nf*z);  uT = zer;
        if direct && src.zs < 0, [uL, uT, uz] = addDirect(uL, uT, uz, m, src, k, phi, w, z, wt); end
    else
        na = m.verticalWavenumber(k, w, 'P');  nb = m.verticalWavenumber(k, w, 'S');
        [uL, uz] = psv(k, na, nb, wt(1)*x.A2.*exp(na*z), wt(2)*x.B2.*exp(nb*z), +1);
        uT = wt(3)*x.C2.*exp(nb*z);
    end
end
end

function [uL, uT, uz] = addDirect(uL, uT, uz, m, src, k, phi, w, z, wt)
% Add the whole-space field e^{-nu|z - zs|} of the source (s = -eps in the table of Sec. 1).
e = sign(z - src.zs);  if e == 0, e = 1; end
d = abs(z - src.zs);
[Pi, Sig, H] = src.directAmplitudes(m, k, phi, w, e);
Pi = wt(1)*Pi;  Sig = wt(2)*Sig;  H = wt(3)*H;
na = m.verticalWavenumber(k, w, 'P');
if m.isFluid
    uL = uL + 1i*k.*Pi.*exp(-na*d);
    uz = uz - e*na.*Pi.*exp(-na*d);
else
    nb = m.verticalWavenumber(k, w, 'S');
    [dL, dz] = psv(k, na, nb, Pi.*exp(-na*d), Sig.*exp(-nb*d), -e);
    uL = uL + dL;  uz = uz + dz;  uT = uT + H.*exp(-nb*d);
end
end

function [uL, uz] = psv(k, na, nb, a, b, s)
% P-SV displacements for phi = a e^{s nu_a z}, psi = b e^{s nu_b z} (table of Sec. 1).
uL = 1i*k.*a - s*nb.*b;
uz = s*na.*a + 1i*k.*b;
end
