function W = interfaceFlux2D(iface, k, omega, res)
%INTERFACEFLUX2D  Horizontal power of an interface-wave field per unit wavefront width.
%   W = lamb.energy.interfaceFlux2D(iface, k, omega, res)
%   res : amplitudes of the mode at the real pole k, fields as iface.unknowns()
%         (.A1 .B1 [.A2 .B2 | .Af]), any common size (e.g. 1-by-Nphi).
%   W   : (omega/2) int Im(s_LL u_L* + s_Lz u_z*) dz over both half-spaces, with the depth
%         integrals done analytically (terms 1/(nu_m + nu_n)). W is a Hermitian quadratic
%         form of the amplitudes: W(c*res) = |c|^2 W(res).
I = halfSpace(iface.lower, k, omega, res.A1, res.B1, -1);
switch iface.kind
    case 'solid-solid'
        I = I + halfSpace(iface.upper, k, omega, res.A2, res.B2, +1);
    case 'solid-fluid'
        f  = iface.upper;
        nf = f.verticalWavenumber(k, omega, 'P');
        I  = I + (-f.rho*omega^2*res.Af).*conj(1i*k*res.Af)/(2*nf);
end
W = 0.5*omega*imag(I);
end

function I = halfSpace(m, k, omega, a, b, s)
% int (s_LL u_L* + s_Lz u_z*) dz for phi = a e^{s nu_a z}, psi = b e^{s nu_b z}
% (s = -1 lower half-space, s = +1 upper). At a real pole every nu is real and positive.
na = m.verticalWavenumber(k, omega, 'P');  nb = m.verticalWavenumber(k, omega, 'S');
G  = 2*k^2 - m.wavenumber(omega, 'S')^2;
ka2 = m.wavenumber(omega, 'P')^2;
uL  = {1i*k*a,                               -s*nb*b};
uz  = {s*na*a,                                1i*k*b};
sLz = {2i*m.mu*k*s*na*a,                     -m.mu*G*b};
sLL = {(-m.lambda*ka2 - 2*m.mu*k^2)*a,       -2i*m.mu*k*s*nb*b};
nu  = [na, nb];
I = 0;
for i = 1:2
    for j = 1:2
        I = I + (sLL{i}.*conj(uL{j}) + sLz{i}.*conj(uz{j}))/(nu(i) + nu(j));
    end
end
end
