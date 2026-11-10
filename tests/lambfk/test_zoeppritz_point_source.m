%TEST_ZOEPPRITZ_POINT_SOURCE  Full-waveform point-source energy vs plane-wave Zoeppritz.
%
%   A deep explosion (H = 100 wavelengths below a welded solid-solid interface) sends a P
%   wavefront that is locally plane at the interface. For incidence angle theta, one receiver
%   is placed on the specular ray of each scattered wave (Rpp, Rps in medium 1; Tpp, Tps in
%   medium 2). The synthetic seismograms (lamb.wavefield.synthetics, DWN + complex frequency)
%   are windowed around each ray arrival. Energy flux I = rho*c*int |v|^2 dt is corrected for
%   geometrical spreading with the exact ray-tube factor of a plane interface and compared with
%   the plane-wave energy coefficients of MaterialClass.ZoeppritzSolid.
%
%   Geometrical spreading (energy conservation in a ray tube; explosion = isotropic source):
%     power into the tube        dP  = I_ref R_ref^2 sin(i1) di1 dphi
%     tube area at the receiver  dA  = X dphi * |dX/di1| cos(i_out) di1
%     =>  E_num = (I_rec / I_ref) * X |dX/di1| cos(i_out) / (R_ref^2 sin(i1)),
%   X = h1 tan(i1) + h2 tan(i_out), sin(i_out) = (c_out/alpha1) sin(i1). In a homogeneous
%   medium this reduces to (I_rec/I_ref)(r/R_ref)^2, the usual r correction.
%   I_ref is measured on the direct P wave of the same source in a whole space at R_ref.
%
%   The explosion radiates only P, so every arrival comes from P incidence (E_Rpp, E_Rpsv,
%   E_Tpp, E_Tpsv). Receivers in medium 1 sit at the source depth, which separates the direct
%   wave from the reflections in time; the direct field is excluded (dwn 'Direct', false).

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));  AddRTMCLib(fullfile(fileparts(mfilename('fullpath')), '..', '..'), false);

assert(exist('MaterialClass', 'class') == 8, 'MaterialClass not found on the path.');

% --- model (shale over sandstone; medium 1 = source side, z > 0) -------------------------
vp1 = 2500; vs1 = 1100; rho1 = 2200;
vp2 = 3000; vs2 = 1700; rho2 = 2300;
m1 = lamb.Medium(rho1, vp1, vs1);  m2 = lamb.Medium(rho2, vp2, vs2);
iface = lamb.Interface(m1, m2);

theta = 30;                                   % incidence angle [deg]
fp    = 25;                                   % Ricker peak frequency [Hz]
lam   = vp1/fp;                               % P wavelength in medium 1
H     = 100*lam;                              % source depth below the interface
t0    = 1.5/fp;                               % wavelet delay
src   = lamb.PointSource.explosion(1e15, H, 'Wavelet', lamb.PointSource.ricker(fp, t0));

% --- receivers on the specular rays -----------------------------------------------------
waves = struct('name', {'Rpp', 'Rps', 'Tpp', 'Tps'}, ...
               'c',    {vp1, vs1, vp2, vs2}, 'rho', {rho1, rho1, rho2, rho2}, ...
               'z',    {H, H, -H, -H});       % medium 1: z = +H (source depth); medium 2: -H
i1 = deg2rad(theta);
for w = 1:numel(waves)
    io = asin(waves(w).c/vp1*sin(i1));        % Snell: outgoing angle
    waves(w).io = io;
    waves(w).X  = H*tan(i1) + H*tan(io);
    dXdi = H/cos(i1)^2 + H/cos(io)^2*(waves(w).c*cos(i1))/(vp1*cos(io));
    waves(w).G  = waves(w).X*abs(dXdi)*cos(io)/sin(i1);    % divided by R_ref^2 below
    waves(w).t  = H/(vp1*cos(i1)) + H/(waves(w).c*cos(io)) + t0;
end
rcv = struct('x', [waves.X].', 'y', zeros(numel(waves), 1), 'z', [waves.z].');

% --- synthetics ----------------------------------------------------------------------------
dt = 0.004;  T = max([waves.t]) + 2;
fprintf('Running DWN synthetics (H = %.0f m = 100 lambda, theta = %g deg)...\n', H, theta);
tic
[t, u] = lamb.wavefield.synthetics(iface, src, rcv, 'T', T, 'dt', dt, 'Direct', false);
Rref = H/2;                                   % reference: direct P in a whole space
rref = struct('x', 0, 'y', 0, 'z', H - Rref);
[~, uref] = lamb.wavefield.synthetics(lamb.Interface(m1, m1), src, rref, 'T', T, 'dt', dt);
fprintf('done in %.0f s\n', toc);

% --- windowed energy flux I = rho c int |v|^2 dt -------------------------------------------
half = 3/fp;                                  % window half-width (Ricker ~ 2/fp long)
flux = @(uu, tc, rho, c) rho*c*trapz(t, sum(gradient(uu.', dt).'.^2 .* ...
                                 (abs(t - tc) < half), 2));
Iref = flux(squeeze(uref(:, :, 1)), Rref/vp1 + t0, rho1, vp1);

% arrivals of the OTHER wave type at each receiver must fall outside the window
other = {'Rps', 'Rpp', 'Tps', 'Tpp'};
for w = 1:numel(waves)
    o = waves(strcmp({waves.name}, other{w}));
    to = rayTime(waves(w).X, H, vp1, H, o.c) + t0;
    assert(abs(to - waves(w).t) > 2*half, '%s window contaminated by %s', waves(w).name, o.name);
end

% --- compare with plane-wave Zoeppritz -------------------------------------------------------
Z = MaterialClass.ZoeppritzSolid(theta, vp1, vs1, rho1, vp2, vs2, rho2);
Ez = 100*[Z.E_Rpp, Z.E_Rpsv, Z.E_Tpp, Z.E_Tpsv];

fprintf('\nP incidence at %g deg, H = 100 lambda, explosion, Ricker %g Hz\n', theta, fp);
fprintf('%-5s %10s %8s %9s | %12s %12s | %9s\n', 'wave', 'offset [m]', 'z [m]', 't [s]', ...
        'numerical %', 'Zoeppritz %', 'rel. diff');
En = zeros(1, numel(waves));
for w = 1:numel(waves)
    Irec  = flux(squeeze(u(:, :, w)), waves(w).t, waves(w).rho, waves(w).c);
    En(w) = 100*Irec/Iref*waves(w).G/Rref^2;
    fprintf('%-5s %10.1f %8.0f %9.3f | %12.4f %12.4f | %8.2f%%\n', waves(w).name, waves(w).X, ...
            waves(w).z, waves(w).t, En(w), Ez(w), 100*abs(En(w) - Ez(w))/Ez(w));
end
fprintf('%-5s %40s | %12.4f %12.4f |\n', 'sum', '', sum(En), sum(Ez));
fprintf('Reflected (P + SV): numerical %.4f %%   Zoeppritz %.4f %%\n', En(1) + En(2), Ez(1) + Ez(2));

tol = 0.005;                                  % 0.5 % relative (finite frequency, 1/(kL) terms)
assert(all(abs(En - Ez)./Ez < tol), 'Numerical energy partition differs from Zoeppritz by > 0.5 %%.');
fprintf('PASS: point-source full-waveform energies match plane-wave Zoeppritz within %g %%.\n', 100*tol);

function t = rayTime(X, h1, v1, h2, v2)
% Travel time of the two-leg ray (h1 at v1, h2 at v2) reaching horizontal offset X.
imax = pi/2 - 1e-9;
if v2 > v1, imax = asin(v1/v2) - 1e-9; end
f  = @(i) h1*tan(i) + h2*tan(asin(v2/v1*sin(i))) - X;
i  = fzero(f, [0, imax]);
io = asin(v2/v1*sin(i));
t  = h1/(v1*cos(i)) + h2/(v2*cos(io));
end
