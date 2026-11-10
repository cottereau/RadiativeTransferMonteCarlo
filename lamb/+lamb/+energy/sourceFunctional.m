function S = sourceFunctional(src, iface, k, phi, omega, x)
%SOURCEFUNCTIONAL  S[u_sc] of the scattered field at the source (eq:Wsrc3).
%   S = lamb.energy.sourceFunctional(src, iface, k, phi, omega, x)
%   force  : S = conj(F) . u(x_s)
%   moment : S = sum_waves v.' * conj(M) * kappa   (= M*_pq d_q u_p)
%   general: the sum of both
%   x holds scattered amplitudes (N-by-Nphi, or residues 1-by-Nphi). Output N-by-Nphi.
k = k(:);  phi = phi(:).';  ik = 1i*k;
l = cat(3, cos(phi), sin(phi), zeros(size(phi)));
t = cat(3, -sin(phi), cos(phi), zeros(size(phi)));
z = cat(3, 0, 0, 1);
if src.zs >= 0                                   % scattered waves e^{-nu z} in the solid
    m = iface.lower;  h = src.zs;
    na = m.verticalWavenumber(k, omega, 'P');  nb = m.verticalWavenumber(k, omega, 'S');
    kP = ik.*l - na.*z;  kS = ik.*l - nb.*z;
    waves = {x.A1.*exp(-na*h), kP,              kP
             x.B1.*exp(-nb*h), nb.*l + ik.*z,   kS};
    if isfield(x, 'C1')
        waves(end+1, :) = {x.C1.*exp(-nb*h), t, kS};
    end
else                                             % reflected fluid wave e^{+nu_f z}
    f = iface.upper;  h = -src.zs;
    nf = f.verticalWavenumber(k, omega, 'P');
    kP = ik.*l + nf.*z;
    waves = {x.Af.*exp(-nf*h), kP, kP};
end
S = 0;
for i = 1:size(waves, 1)
    [a, p, kap] = waves{i, :};
    if ~strcmp(src.kind, 'moment')
        S = S + a.*sum(p.*reshape(conj(src.F), 1, 1, 3), 3);
    end
    if ~strcmp(src.kind, 'force')
        for pp = 1:3
            for qq = 1:3
                S = S + conj(src.M(pp, qq))*a.*p(:, :, pp).*kap(:, :, qq);
            end
        end
    end
end
end
