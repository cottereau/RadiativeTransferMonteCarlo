function Wi = interfaceWaves(iface, src, poles, omega, opts)
%INTERFACEWAVES  Power trapped in each real pole, computed two independent ways.
%   Wi = lamb.energy.interfaceWaves(iface, src, poles, omega, opts)
%
%   (a) far-field flux (eq:WIW3):  W = k_p/(2 pi) int_0^{2pi} W2D+(phi) dphi, where W2D+ is
%       the horizontal flux (omega/2) int Im(s_LL u_L* + s_Lz u_z*) dz of the residue fields,
%       with the depth integrals done analytically (terms 1/(nu_m + nu_n)).
%   (b) source side (last term of eq:Wsrc3): Wsrc = omega k_p/(8 pi) int Re S[Res u_sc] dphi.
%   Output struct array: .type .k .c .W (a) .Wsrc (b)
[phi, wphi] = lamb.quad.azimuthNodes(opts.NPhi);
Wi = struct('type', {}, 'k', {}, 'c', {}, 'W', {}, 'Wsrc', {});
for p = poles
    res = lamb.poleResidues(iface, src, p, phi, omega);
    W2D = lamb.energy.interfaceFlux2D(iface, p.k, omega, res);
    S   = lamb.energy.sourceFunctional(src, iface, p.k, phi, omega, res);
    Wi(end+1) = struct('type', p.type, 'k', p.k, 'c', p.c, ...
        'W',    p.k/(2*pi)*sum(wphi.*W2D), ...
        'Wsrc', omega*p.k/(8*pi)*sum(wphi.*real(S)));             %#ok<AGROW>
end
end
