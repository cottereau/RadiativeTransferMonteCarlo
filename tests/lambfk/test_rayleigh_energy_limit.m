% References:
%   - Miller, G.F. & Pursey, H. (1955), On the partition of energy between elastic waves in
%     a semi-infinite solid, Proc. R. Soc. Lond. A 233, 55-69. doi:10.1098/rspa.1955.0245

%TEST_RAYLEIGH_ENERGY_LIMIT  Solid-fluid interface -> free surface: energy partition check.
%
%   Poisson solid (lambda = mu, alpha = sqrt(3) beta) under a vanishing fluid
%   (rho_f = alpha_f = 1e-6), vertical point force at the interface (h = 0).
%   Expected (Miller & Pursey 1955): Rayleigh ~67.4 %, S ~25.8 %, P ~6.9 %.
%
%   Physics note. For alpha_f < c_R the Rayleigh pole of a solid-fluid interface is LEAKY
%   (it radiates into the fluid): k_L = k_R + 1i*b with b proportional to rho_f*alpha_f
%   (eq:dk). Here b ~ 1e-13, i.e. a leakage length 1/(2b) ~ 1e12 wavelengths: the wave is
%   trapped for every practical range. Its energy appears in the fluid-P window as a
%   Lorentzian peak of half-width b, which the quadrature resolves with graded nodes;
%   EP.pct.trapped = real interface poles + quasi-trapped leaky poles.
%   For comparison, the same limit is also run with a FAST fluid (alpha_f = 1.5 > c_R),
%   where the trapped pole is a true (real) Scholte pole that tends to the Rayleigh pole.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));  AddRTMCLib(fullfile(fileparts(mfilename('fullpath')), '..', '..'), false);

% --- 1. media ------------------------------------------------------------------------
beta = 1;  rho = 1;
solid = lamb.Medium(rho, sqrt(3)*beta, beta);                % lambda = mu
fprintf('Poisson solid: lambda/mu = %.6f\n', solid.lambda/solid.mu);

% --- 2. fluids approaching vacuum --------------------------------------------------------
fluids = {lamb.Medium(1e-6, 1e-6), lamb.Medium(1e-6, 1.5)};
labels = {'rho_f = 1e-6, alpha_f = 1e-6 (slow: leaky Rayleigh)', ...
          'rho_f = 1e-6, alpha_f = 1.5  (fast: real Scholte)'};

% --- 3. point source at the boundary (h = 0), vertical force ----------------------------
src   = lamb.PointSource('force', [0; 0; 1], 0);
omega = 2*pi;                                                % any omega: h = 0 has no scale

% reference: free surface
EPfs = lamb.energyPartition(lamb.Interface(solid), src, omega);
ref  = [EPfs.pct.body.P1, EPfs.pct.body.SV1 + EPfs.pct.body.SH1, EPfs.pct.ifTotal];
fprintf('\n%-52s  P [%%]   S [%%]  Fluid [%%]  Trapped [%%]  balance err\n', 'Model');
fprintf('%-52s  %6.2f  %6.2f  %8s  %10.2f   %8.1e\n', 'free surface (reference)', ...
        ref(1), ref(2), '-', ref(3), EPfs.check.relErr);

% --- 4. energy partition for the fluid-loaded cases -------------------------------------
tol = 0.05;                                                  % percentage points
mp  = [6.9 25.8 67.4];                                       % Miller & Pursey (1955)
for i = 1:numel(fluids)
    EP = lamb.energyPartition(lamb.Interface(solid, fluids{i}), src, omega);
    P  = EP.pct.body.P1;
    S  = EP.pct.body.SV1 + EP.pct.body.SH1;
    Fp = EP.pct.body.Pf - sum(EP.pct.leaky);                 % fluid P outside the leaky peak
    R  = EP.pct.trapped;
    fprintf('%-52s  %6.2f  %6.2f  %8.4f  %10.2f   %8.1e\n', labels{i}, P, S, Fp, R, EP.check.relErr);
    for L = EP.W.leaky
        fprintf('    %s pole: k = %.10f %+.3ei, c = %.6f beta, leakage length = %.2e wavelengths\n', ...
                L.type, real(L.k), imag(L.k), omega/real(L.k), L.leakLength*real(L.k)/(2*pi));
    end
    for p = EP.W.interface
        fprintf('    %s pole: c = %.6f beta\n', p.type, p.c);
    end

    % --- 5. checks -----------------------------------------------------------------------
    assert(abs(R - ref(3)) < tol,         'Trapped energy %.3f %% differs from Rayleigh %.3f %%', R, ref(3));
    assert(abs(P + S - (ref(1) + ref(2))) < tol, 'Body waves %.3f %% differ from %.3f %%', P + S, ref(1) + ref(2));
    assert(abs(R - mp(3)) < 0.1 && abs(P + S - (mp(1) + mp(2))) < 0.1, 'Not the Miller-Pursey partition');
    assert(EP.check.relErr < 1e-6,        'Energy balance error %.2e', EP.check.relErr);
end
fprintf('\nPASS: interface wave -> Rayleigh %.2f %% (Miller & Pursey 67.4 %%), body waves %.2f %% (~33 %%)\n', ...
        ref(3), ref(1) + ref(2));
