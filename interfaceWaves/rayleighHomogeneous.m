function R = rayleighHomogeneous(geometry,source,material,observation,opts)
%RAYLEIGHHOMOGENEOUS  Rayleigh-wave channel for a homogeneous cylinder.
%
%   R = rayleighHomogeneous(geometry,source,material,observation,opts)
%   propagates the Rayleigh waves launched by a source sitting on a free
%   surface of a homogeneous cylinder, and returns their energy on exactly
%   the grid radiativeTransfer uses, so R.energyDensity can be added to or
%   overlaid on obs.energyDensity.
%
%   WHY THIS IS A SEPARATE FUNCTION. In a homogeneous medium the Rayleigh
%   and body-wave populations are decoupled. A flat free surface cannot
%   turn a body wave into a Rayleigh wave: the Rayleigh horizontal slowness
%   1/cR exceeds 1/vS, so no body wave can phase match it. The coupling of
%   Xu et al. (2021) is carried entirely by the scatterers and vanishes
%   with them. What is left is the source itself, plus the edges.
%
%   The source partition defaults to a normal point load on a free surface,
%     Miller, G.F. & Pursey, H. (1955). On the partition of energy between
%     elastic waves in a semi-infinite solid. Proc. R. Soc. Lond. A 233,
%     55-69. doi:10.1098/rspa.1955.0245
%   which gives 67.4 per cent Rayleigh, 25.8 per cent shear and 6.9 per
%   cent dilatational for a Poisson solid (nu = 0.25). Override it through
%   opts.partition for a different Poisson ratio or a real transducer.
%
%   At each 90 degree edge the wave splits into a part that travels round
%   onto the neighbouring surface, a part reflected back and a part
%   radiated into the bulk. The defaults are order-of-magnitude values for
%   a right-angle edge in a metal (see e.g. Viktorov, I.A., 1967, Rayleigh
%   and Lamb Waves, Plenum Press, ch. 1); they are not settled to better
%   than tens of per cent, so calibrate opts.edge against your own data.
%   The radiated share is removed from the Rayleigh population and tallied
%   in R.radiated rather than fed back into the body-wave field.
%
%   opts (all optional)
%     numberParticles  default 2e5. The Rayleigh front is coherent, so it
%                      needs far fewer particles than the diffuse body field
%     partition        [fP fS fR] source energy split, default [.069 .258 .674]
%     edge             [transmitted reflected] energy fractions at an edge,
%                      default [0.30 0.12]; the rest radiates into the bulk
%     direction        'uniform' or 'outgoing'; defaults to source.direction
%     substeps         sub-steps per observation time step, default 4
%
%   See also RAYLEIGHMODEL, RADIATIVETRANSFER, TESTRAYLEIGHHOMOGENEOUS.

if nargin < 5, opts = struct; end

% Formal 90-degree wedge coefficients based on material Poisson's ratio
if isfield(material, 'vp') && isfield(material, 'vs') && ~isempty(material.vp) && ~isempty(material.vs)
    c_wedge = rayleighWedgeCoefficients([], material.vp, material.vs);
    edge_default = [c_wedge.TR, c_wedge.RR];
else
    edge_default = [0.38 0.14]; % Poisson solid default
end

def = struct('numberParticles',2e5,'partition',[0.069 0.258 0.674], ...
             'edge',edge_default,'substeps',4,'direction','');
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts,f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
end
if isempty(opts.direction)
    if isfield(source,'direction') && ~isempty(source.direction)
        opts.direction = source.direction;
    else
        opts.direction = 'uniform';
    end
end

%% 1. geometry: a cylinder closed by two planes normal to z
ztop = []; zbot = []; rmax = [];
zs = [];
for i = 1:numel(geometry.bnd)
    b = geometry.bnd(i);
    if b.dir == 3,     zs(end+1) = b.val; %#ok<AGROW>
    elseif b.dir == 4, rmax = b.val;
    end
end
if numel(zs) ~= 2 || isempty(rmax)
    error(['rayleighHomogeneous expects a cylinder: two boundaries with ' ...
           'dir = 3 and one with dir = 4.']);
end
ztop = max(zs);  zbot = min(zs);
if ~strcmp(geometry.frame,'cylindrical')
    error('rayleighHomogeneous expects geometry.frame = ''cylindrical''.');
end

%% 2. Rayleigh velocity and depth profile
M  = RayleighModel(material.vp,material.vs,material.Frequency,0,0,1,ztop-zbot);
cR = M.cR;
lambdaR = cR/material.Frequency;

% G(z) = int_0^z (r1^2+r2^2) dz' is already tabulated (with n = 1). It is
% normalised by its own total rather than by M.Inorm so that the bin shares
% sum to one exactly: the two differ by the quadrature error of the table.
Gtab = M.G1 + M.G2;
G = @(z) interp1(M.zgrid,Gtab,min(max(z,0),ztop-zbot),'linear');
W = @(d0,d1) max(0,(G(d1)-G(d0))/Gtab(end)); % energy share between two depths

%% 3. observation grid, matching initializeObservation / radiativeTransfer
binX = observation.x(:).';   xc = (binX(1:end-1)+binX(2:end))/2;  Nx = numel(xc);
binZ = observation.z(:).';                                        Nz = numel(binZ)-1;
tobs = [0 setdiff(observation.time,0)];                           Nt = numel(tobs);
dxv  = xc.*diff(binX);          % cylindrical volume weight, r dr
dzv  = diff(binZ);
dyv  = 2*pi;                    % observation.y spans the full azimuth

% energy shares over the bins, one vector per surface
wTop = W(ztop-binZ(2:end), ztop-binZ(1:end-1));        % depth measured downwards
wBot = W(binZ(1:end-1)-zbot, binZ(2:end)-zbot);        % depth measured upwards
wLat = W(rmax-binX(2:end), rmax-binX(1:end-1));        % depth measured inwards
if lambdaR > 0.25*(ztop-zbot)
    warning(['lambda_R = %.3g is not small against the height %.3g: the two ' ...
             'faces couple and the modes are Lamb-like, not Rayleigh.'], ...
             lambdaR,ztop-zbot);
end

%% 4. source on the top free surface
N  = round(opts.numberParticles);
fR = opts.partition(3);
if abs(source.position(3)-ztop) > lambdaR
    warning(['the source is %.3g from the free surface, more than one ' ...
             'Rayleigh wavelength: direct Rayleigh excitation is not ' ...
             'physical there.'],abs(source.position(3)-ztop));
end

r0 = zeros(N,1);
if isfield(source,'radial') && ~isempty(source.radial)
    % footprint on the surface: pdf(r) ~ amp(r) * r  (a disc, not a ball)
    rg  = source.radial.GridVectors{1}(:);
    pdf = max(0,source.radial(rg)).*rg;
    cdf = cumtrapz(rg,pdf);  cdf = cdf/cdf(end);
    keep = [true; diff(cdf) > 0];
    r0 = interp1(cdf(keep),rg(keep),rand(N,1),'linear','extrap');
end
th0 = 2*pi*rand(N,1);
P.surf = ones(N,1);                       % 1 top face, 2 bottom face, 3 lateral
P.u = r0.*cos(th0);
P.v = r0.*sin(th0);
if strcmpi(opts.direction,'outgoing')
    ph = th0;
else
    ph = 2*pi*rand(N,1);
end
P.du = cos(ph);  P.dv = sin(ph);
P.alive = true(N,1);

%% 5. propagate and record
nsub = max(1,round(opts.substeps));

E    = zeros(Nx,Nz,Nt);
rad  = zeros(Nt,1);
lost = 0;
Tr   = opts.edge(1);  Rf = opts.edge(2);
if Tr+Rf > 1, error('opts.edge must sum to at most one.'); end

for it = 1:Nt
    if it > 1
        % Advance by this interval, not by mean(diff(tobs)): observation.time
        % is not required to be uniform, and a mean step silently misplaces
        % every particle relative to the time it is recorded at.
        step = cR*(tobs(it)-tobs(it-1))/nsub;
        for is = 1:nsub
            P.u = P.u + step*P.du;
            P.v = P.v + step*P.dv;
            [P,dl] = edges(P,rmax,ztop,zbot,Tr,Rf);
            lost = lost + dl;
        end
    end
    rad(it) = lost/N;

    a = P.alive;
    for sf = 1:3
        k = a & P.surf == sf;
        if ~any(k), continue; end
        if sf < 3
            ir = discretize(snap(hypot(P.u(k),P.v(k)),binX),binX);
            ir = ir(~isnan(ir));
            if isempty(ir), continue; end
            c  = accumarray(ir,1,[Nx 1]);
            if sf == 1, E(:,:,it) = E(:,:,it) + c*wTop;
            else,       E(:,:,it) = E(:,:,it) + c*wBot;
            end
        else
            iz = discretize(snap(P.v(k),binZ),binZ);
            iz = iz(~isnan(iz));
            if isempty(iz), continue; end
            c  = accumarray(iz,1,[Nz 1]);
            E(:,:,it) = E(:,:,it) + wLat.'*c.';
        end
    end
end

%% 6. normalise exactly like radiativeTransfer
E = E*(fR/N);                                  % energy per source unit
R.energy        = squeeze(sum(sum(E,1),2));    % Rayleigh energy vs time
R.energyDensity = E./(dxv.'*dzv)/dyv;
R.radiated      = rad*fR;                      % energy given to the bulk at edges
R.t = tobs;  R.x = xc;  R.z = (binZ(1:end-1)+binZ(2:end))/2;
R.cR = cR;  R.lambdaR = lambdaR;  R.M = M;
R.partition = opts.partition;  R.edge = [Tr Rf 1-Tr-Rf];
R.rmax = rmax;  R.ztop = ztop;  R.zbot = zbot;
R.depthWeights = struct('top',wTop,'bottom',wBot,'lateral',wLat);

end

% -------------------------------------------------------------------------
function v = snap(v,edges)
% a particle pinned to a boundary lands one ulp outside the outer bin edge,
% which discretize reports as out of range. Pull those back in; anything
% genuinely outside the observation grid is left alone.
lo = edges(1);  hi = edges(end);
v(v < lo & v > lo-eps(lo)*8) = lo;
v(v > hi & v < hi+eps(hi)*8) = hi;
end

% -------------------------------------------------------------------------
function [P,lost] = edges(P,rmax,ztop,zbot,Tr,Rf)
% 90 degree edges: travel round onto the next surface, reflect, or radiate
lost = 0;

% --- a face particle reaching the rim ---
for sf = 1:2
    r = hypot(P.u,P.v);
    k = find(P.alive & P.surf == sf & r > rmax);
    if isempty(k), continue; end
    er = [P.u(k) P.v(k)]./r(k);
    et = [-er(:,2) er(:,1)];
    dr = P.du(k).*er(:,1) + P.dv(k).*er(:,2);     % outward component
    dt = P.du(k).*et(:,1) + P.dv(k).*et(:,2);     % along the edge
    P.u(k) = er(:,1)*rmax;  P.v(k) = er(:,2)*rmax;

    q = rand(numel(k),1);
    go = q < Tr;  bk = ~go & q < Tr+Rf;  kill = ~go & ~bk;

    if any(go)                                    % round onto the lateral face
        g = k(go);
        P.surf(g) = 3;
        P.u(g) = rmax*atan2(er(go,2),er(go,1));   % arc length
        if sf == 1, P.v(g) = ztop; P.dv(g) = -dr(go);
        else,       P.v(g) = zbot; P.dv(g) =  dr(go);
        end
        P.du(g) = dt(go);
    end
    if any(bk)                                    % reflect back onto the face
        g = k(bk);
        P.du(g) = -dr(bk).*er(bk,1) + dt(bk).*et(bk,1);
        P.dv(g) = -dr(bk).*er(bk,2) + dt(bk).*et(bk,2);
    end
    if any(kill)
        P.alive(k(kill)) = false;
        lost = lost + nnz(kill);
    end
end

% --- a lateral particle reaching either rim ---
for e = 1:2
    if e == 1, k = find(P.alive & P.surf == 3 & P.v > ztop);
    else,      k = find(P.alive & P.surf == 3 & P.v < zbot);
    end
    if isempty(k), continue; end
    dz = P.dv(k);  ds = P.du(k);
    if e == 1, P.v(k) = ztop; else, P.v(k) = zbot; end

    q = rand(numel(k),1);
    go = q < Tr;  bk = ~go & q < Tr+Rf;  kill = ~go & ~bk;

    if any(go)                                    % round onto a face
        g = k(go);
        th = P.u(g)/rmax;
        er = [cos(th) sin(th)];  et = [-sin(th) cos(th)];
        P.surf(g) = 1 + (e == 2);
        P.u(g) = rmax*er(:,1);  P.v(g) = rmax*er(:,2);
        P.du(g) = -abs(dz(go)).*er(:,1) + ds(go).*et(:,1);   % inward
        P.dv(g) = -abs(dz(go)).*er(:,2) + ds(go).*et(:,2);
    end
    if any(bk), P.dv(k(bk)) = -dz(bk); end
    if any(kill)
        P.alive(k(kill)) = false;
        lost = lost + nnz(kill);
    end
end
end
