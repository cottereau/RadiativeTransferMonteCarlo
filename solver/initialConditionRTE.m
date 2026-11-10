function IC = initialConditionRTE(ampTabFile,vp,t0,Nr)
%INITIALCONDITIONRTE  RTE initial condition for an isotropic moment source.
%
%   IC = initialConditionRTE(ampTabFile,vp,t0,Nr) turns the SEM source time
%   function into the radial profile the Monte Carlo code expects.
%
%   An isotropic moment tensor M_ij = M0 g(t) delta_ij delta(x) radiates P
%   only, with no S at all in a full space. Its far field is
%
%       u_r = M0 gdot(t-r/alpha) / (4 pi rho alpha^3 r)
%       E(r,t) = rho |udot_r|^2 = M0^2 gddot(t-r/alpha)^2 /(16 pi^2 rho alpha^6 r^2)
%
%   so at a reference time t0 after the source has stopped, the number of
%   particles per unit radius goes as gddot(t0-r/alpha)^2: the r^2 volume
%   weight cancels the 1/r^2 of the spherical spreading exactly.
%
%   CONVENTION. initializeParticle feeds source.radial to inverseCDF, which
%   applies the r^(d-1) weight itself (initializeParticle.m:184-188). So
%   source.radial must be the ENERGY DENSITY, carrying the 1/r^2, and not
%   the radial pdf. IC.amp is built that way; IC.pdf is the resulting
%   particle density per unit radius, for plotting and for checks.
%
%   Two things that are easy to get wrong:
%     - the map is RETARDED, r = alpha*(t0-tau): what left the source first
%       sits at the LARGEST radius, so the profile is time reversed;
%     - the RTE clock then starts at t0, so traces must be shifted by t0,
%       not by the source duration.
%
%   ampTabFile is the tabulated moment function M(t) (input.spec: func =
%   file, time_file). If a run instead tabulates the moment RATE,
%   differentiate once less.
%
%   See also INITIALIZEPARTICLE, RAYLEIGHHOMOGENEOUS.

if nargin < 2 || isempty(vp), vp = 4814.4; end
if nargin < 3 || isempty(t0), t0 = 6e-6;   end
if nargin < 4 || isempty(Nr), Nr = 400;    end

A = readmatrix(ampTabFile);
t = A(:,1);  M = A(:,2);
dt = mean(diff(t));

% Mddot: the far-field particle velocity waveform
Md  = gradient(M,dt);
Mdd = gradient(Md,dt);

% support of the source
act = find(abs(Mdd) > 1e-6*max(abs(Mdd)));
T   = t(act(end));
if t0 < T
    warning(['t0 = %.2f us is inside the source (it stops at %.2f us): the ' ...
             'shell is still attached to the source and is not a free wave.'], ...
             1e6*t0,1e6*T);
end

r   = linspace(0,vp*t0,Nr).';
tau = t0 - r/vp;                              % retarded emission time
g2  = interp1(t,Mdd,tau,'linear',0).^2;
g2(tau < 0 | tau > t(end)) = 0;

amp = g2./max(r,eps).^2;                      % energy density, for source.radial
amp(1) = 0;                                   % the shell never reaches r = 0

IC.r = r;  IC.amp = amp;  IC.pdf = g2;
IC.t0 = t0;  IC.vp = vp;  IC.T = T;
IC.rInner = vp*max(0,t0-T);  IC.rOuter = vp*t0;
IC.stf = struct('t',t,'M',M,'Mdot',Md,'Mddot',Mdd);
[~,ip] = max(g2);  IC.rPeak = r(ip);
end
