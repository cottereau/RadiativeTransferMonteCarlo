function rt = rtCoefficients(iface, incident, p, omega)
%RTCOEFFICIENTS  Plane-wave energy reflection/transmission coefficients vs horizontal slowness.
%   rt = lamb.energy.rtCoefficients(iface, 'P1' | 'SV1' | 'SH1' | 'Pf', p, omega)
%
%   p  : horizontal slowness [s/m] (k = omega*p), any shape (returned as a column)
%   rt : .p, one field per outgoing wave (e.g. .P1 .SV1 .SH1 .P2 .SV2 .SH2 .Pf) holding
%        e_m = E_m |x_m|^2 / (E_inc |a_inc|^2), and .sum (= 1 for propagating incidence).
%   'P1','SV1','SH1' are incident from the lower solid; 'Pf' from the fluid.
k = omega*p(:);  N = numel(k);  o = ones(N, 1);  z = zeros(N, 1);
switch incident
    case 'P1',  amp = struct('P0', o, 'S0', z, 'C0', z); type = 'P';
    case 'SV1', amp = struct('P0', z, 'S0', o, 'C0', z); type = 'SV';
    case 'SH1', amp = struct('P0', z, 'S0', z, 'C0', o); type = 'SH';
    case 'Pf',  amp = struct('Pf0', o);                   type = 'P';
    otherwise,  error('lamb:rtCoefficients:incident', 'Unknown incident wave %s.', incident);
end
srcInSolid = ~strcmp(incident, 'Pf');
if srcInSolid, medInc = iface.lower; else, medInc = iface.upper; end
Einc = lamb.energy.fluxDensity(medInc, k, omega, type);
x = iface.solve(amp, k, omega);
rt.p = p(:);  rt.sum = 0;
for w = lamb.energy.outgoingWaves(iface, srcInSolid)
    rt.(w.name) = lamb.energy.fluxDensity(w.med, k, omega, w.type).*abs(x.(w.field)).^2./Einc;
    rt.sum = rt.sum + rt.(w.name);
end
end
