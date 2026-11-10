function amp = incidentAmplitudes(src, iface, k, phi, omega, away)
%INCIDENTAMPLITUDES  Source field decomposed into plane waves at the interface (eq:PSH).
%   amp = lamb.incidentAmplitudes(src, iface, k, phi, omega)
%   amp = lamb.incidentAmplitudes(src, iface, k, phi, omega, true)   % also 'away' fields
%
%   k N-by-1, phi 1-by-Nphi, omega scalar (real or complex). All fields are N-by-Nphi.
%   Source in the solid (zs = h >= 0):
%     .P0 .S0 .C0      toward the interface, proportional to e^{+nu z} (enter M x = S)
%     .Pd0 .Sd0 .Cd0   away from the interface, proportional to e^{-nu z}   (if away)
%   Source in the fluid (zs = -h < 0):
%     .Pf0             toward the interface, proportional to e^{-nu_f z}
%     .Uf0             away (upward), proportional to e^{+nu_f z}           (if away)
%   The 'away' fields contain e^{+nu h}: use them only inside the propagating windows.
arguments
    src   (1,1) lamb.PointSource
    iface (1,1) lamb.Interface
    k     (:,1) double
    phi   (1,:) double
    omega (1,1) double
    away  (1,1) logical = false
end
if src.zs >= 0
    med = iface.lower;  h = src.zs;
    na = med.verticalWavenumber(k, omega, 'P');
    nb = med.verticalWavenumber(k, omega, 'S');
    [Pi, Sig, H] = src.directAmplitudes(med, k, phi, omega, -1);
    amp.P0 = Pi.*exp(-na*h);  amp.S0 = Sig.*exp(-nb*h);  amp.C0 = H.*exp(-nb*h);
    if away
        [Pi, Sig, H] = src.directAmplitudes(med, k, phi, omega, +1);
        amp.Pd0 = Pi.*exp(na*h);  amp.Sd0 = Sig.*exp(nb*h);  amp.Cd0 = H.*exp(nb*h);
    end
else
    if ~strcmp(iface.kind, 'solid-fluid')
        error('lamb:incidentAmplitudes:side', 'zs < 0 requires a fluid above the interface.');
    end
    med = iface.upper;  h = -src.zs;
    nf = med.verticalWavenumber(k, omega, 'P');
    amp.Pf0 = src.directAmplitudes(med, k, phi, omega, +1).*exp(-nf*h);
    if away
        amp.Uf0 = src.directAmplitudes(med, k, phi, omega, -1).*exp(nf*h);
    end
end
end
