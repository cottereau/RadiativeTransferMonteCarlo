function K = slabKernels(T, inc, p, omega, z, het, opts)
%SLABKERNELS  Scattering cross sections of a point heterogeneity near an interface (K2, K3).
%   K = lamb.rte.slabKernels(T, inc, p, omega, z, het)
%   K = lamb.rte.slabKernels(T, inc, p, omega, z, het, 'Mode', 'incident', 'Tflip', Tf)
%
%   T     radiation tables of the interface (lamb.rte.radiationTables(T.iface))
%   inc   incident wave (lamb.rte.scatteringState): 'P1' 'SV1' 'SH1' 'Pf' 'P2' 'SV2' 'SH2'
%         or 'IWj'; p its horizontal slowness (ignored for 'IWj')
%   z     scatterer depths (vector); het  heterogeneity struct (lamb.rte.equivalentSource)
%   Mode  'state' (default, DWBA) or 'incident' (incident wave only at the scatterer)
%   Tflip tables of the flipped interface; needed for scatterers in an UPPER SOLID (z < 0),
%         whose outgoing radiation is computed in the mirrored problem
%
%   The scatterer re-radiates as the equivalent point source s(z) with the half-space
%   Green's function, so each outgoing channel c gets W_c = s' Q_c(z) s (lamb.rte.channelQ).
%   K.sigma.(c)  1-by-Nz cross sections W_c/I (body channels, in the ORIGINAL naming:
%                P1 SV1 SH1 | P2 SV2 SH2 | Pf); for 'IWj' incidence these are rates per
%                unit path length of the mode per unit scatterer density
%   K.iw         nIW-by-Nz cross sections into the interface waves
%   K.leaky      nLeaky-by-Nz, the part of K.sigma.Pf inside the leaky peaks
%   K.total      1-by-Nz sum over body channels and interface waves
%   K.z, K.I, K.s (9-by-Nz equivalent sources)
arguments
    T struct
    inc char
    p double
    omega (1,1) double
    z double
    het struct
    opts.Mode char = 'state'
    opts.Tflip = []
end
z = z(:).';
S = lamb.rte.scatteringState(T.iface, inc, p, omega, z, 'Mode', opts.Mode);
s = lamb.rte.equivalentSource(S, het, omega);
kb1 = omega/T.iface.lower.beta;
K.z = z;  K.I = S.I;  K.s = s;  K.iw = [];  K.leaky = [];
ss = strcmp(T.iface.kind, 'solid-solid');
for iz = 1:numel(z)
    if ss && z(iz) < 0                                   % scatterer in the upper solid
        if isempty(opts.Tflip)
            error('lamb:slabKernels:flip', 'Scatterers in an upper solid need opts.Tflip.');
        end
        sm = s(:, iz).*[1 1 -1 1 1 1 1 -1 -1].';         % mirror z: F_z, M_xz, M_yz change sign
        [~, W] = lamb.rte.channelQ(opts.Tflip, -kb1*z(iz)*T.iface.lower.beta/opts.Tflip.iface.lower.beta, ...
                                   omega, sm);
        W.body = swapNames(W.body);
    else
        [~, W] = lamb.rte.channelQ(T, kb1*z(iz), omega, s(:, iz));
    end
    for f = fieldnames(W.body).'
        K.sigma.(f{1})(iz) = W.body.(f{1})/S.I;
    end
    K.iw(:, iz) = W.iw(:)/S.I;
    K.leaky(:, iz) = W.leaky(:)/S.I;
end
K.total = sum(cell2mat(struct2cell(K.sigma)), 1) + sum(K.iw, 1);
end

function B = swapNames(A)
% channels of the flipped problem: its medium 1 is our medium 2
map = struct('P1', 'P2', 'SV1', 'SV2', 'SH1', 'SH2', 'P2', 'P1', 'SV2', 'SV1', 'SH2', 'SH1');
B = struct();
for f = fieldnames(A).'
    B.(map.(f{1})) = A.(f{1});
end
end
