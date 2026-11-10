function E = fluxDensity(med, k, omega, wave)
%FLUXDENSITY  Vertical power density of one plane wave per unit |amplitude|^2.
%   E = lamb.energy.fluxDensity(med, k, omega, 'P' | 'SV' | 'SH')
%   P, SV (potential amplitude)   : E = 0.5*rho*omega^3*Re(1i*nu)     (eq:Em)
%   SH    (displacement amplitude): E = 0.5*omega*mu*Re(1i*nu_beta)
%   Zero outside the propagating window |k| < k_m (real omega).
switch upper(wave)
    case 'P',  E = 0.5*med.rho*omega^3*real(1i*med.verticalWavenumber(k, omega, 'P'));
    case 'SV', E = 0.5*med.rho*omega^3*real(1i*med.verticalWavenumber(k, omega, 'S'));
    case 'SH', E = 0.5*omega*med.mu   *real(1i*med.verticalWavenumber(k, omega, 'S'));
    otherwise, error('lamb:fluxDensity:wave', 'wave must be P, SV or SH.');
end
end
