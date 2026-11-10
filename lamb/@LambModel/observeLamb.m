function E = observeLamb(M, P)
%OBSERVELAMB  Energy of Lamb-model particles per channel.
%   E = observeLamb(M, P)
%   E.(name) for the body channels of M.chan (P1, SV1, SH1 | P2, SV2, SH2 | Pf), .IW (vector,
%   one per real interface-wave pole), .LIW (vector, one per leaky pole), .total, and .pct (the
%   same fields in % of .total). An S particle is split between SV and SH by aSV^2 and aSH^2 of
%   its polarisation with respect to the interface normal (lambConventionMap).
%
%   See also LAMBRTEMODEL, INITIALIZEPARTICLESLAMB, PROPAGATEPARTICLELAMB.
a = P.alive;
L = lambConventionMap(P.dir, P.perp, [0 0 1]);
for ch = M.chan
    m = a & P.side == ch.side & P.mode == ch.mode;
    switch ch.name(1:2)
        case 'SV', E.(ch.name) = sum(P.E(m).*L.aSV(m).^2);
        case 'SH', E.(ch.name) = sum(P.E(m).*L.aSH(m).^2);
        otherwise, E.(ch.name) = sum(P.E(m));
    end
end
E.IW  = arrayfun(@(j) sum(P.E(a & P.mode == 10 + j)), 1:numel(M.iw));
E.LIW = arrayfun(@(j) sum(P.E(a & P.mode == 20 + j)), 1:numel(M.leaky));
E.total = sum(P.E(a));
E.pct = structfun(@(v) 100*v/E.total, rmfield(E, 'total'), 'UniformOutput', false);
end
