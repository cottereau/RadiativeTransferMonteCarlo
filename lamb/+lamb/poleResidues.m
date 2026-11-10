function res = poleResidues(iface, src, pole, phi, omega)
%POLERESIDUES  Residues of the scattered P-SV amplitudes at a real pole (eq:res).
%   res = lamb.poleResidues(iface, src, pole, phi, omega)
%
%   Res x = r (l.' S(k_p)) / (l.' M'(k_p) r), where r and l are the right and left null
%   vectors of M(k_p) and M' its derivative (lamb.poleNullVectors: Richardson-extrapolated
%   centred difference). M depends on k only, so one SVD serves all azimuths. Output fields are 1-by-Nphi (.A1 .B1 [.A2 .B2 | .Af]).
%   No SH field: two half-spaces have no SH pole.
kp = pole.k;
[r, lT, den] = lamb.poleNullVectors(iface, kp, omega);
amp = lamb.incidentAmplitudes(src, iface, kp, phi, omega);
S  = iface.rhs(kp, omega, amp);                    % n-by-Nphi(-by-1)
R  = r*(lT*S)/den;
names = iface.unknowns();
for i = 1:numel(names)
    res.(names{i}) = R(i, :);
end
end
