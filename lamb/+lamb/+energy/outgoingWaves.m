function w = outgoingWaves(iface, srcInSolid)
%OUTGOINGWAVES  Body waves leaving the interface region (table in sec:energy + SH rows).
%   w = lamb.energy.outgoingWaves(iface, srcInSolid)
%   Struct array with fields
%     .name    'P1' 'SV1' 'SH1' 'P2' 'SV2' 'SH2' 'Pf'
%     .med     lamb.Medium the wave travels in
%     .type    'P' | 'SV' | 'SH'  (for lamb.energy.fluxDensity)
%     .field   scattered amplitude in x (from Interface.solve)
%     .away    direct-field amplitude travelling the same way ('' if none)
%     .toward  direct-field amplitude travelling the other way (for the whole-space power)
L = iface.lower;  U = iface.upper;
if srcInSolid
    w = [wave('P1', L, 'P', 'A1', 'Pd0', 'P0'), wave('SV1', L, 'SV', 'B1', 'Sd0', 'S0'), ...
         wave('SH1', L, 'SH', 'C1', 'Cd0', 'C0')];
    switch iface.kind
        case 'solid-solid'
            w = [w, wave('P2', U, 'P', 'A2', '', ''), wave('SV2', U, 'SV', 'B2', '', ''), ...
                    wave('SH2', U, 'SH', 'C2', '', '')];
        case 'solid-fluid'
            w = [w, wave('Pf', U, 'P', 'Af', '', '')];
    end
else
    w = [wave('Pf', U, 'P', 'Af', 'Uf0', 'Pf0'), wave('P1', L, 'P', 'A1', '', ''), ...
         wave('SV1', L, 'SV', 'B1', '', '')];
end
end

function s = wave(name, med, type, field, away, toward)
s = struct('name', name, 'med', med, 'type', type, 'field', field, ...
           'away', away, 'toward', toward);
end
