function P = propagateParticleMultiMaterial(matArray,geometry,P,T)
%PROPAGATEPARTICLEMULTIMATERIAL  RTE random walk across several materials.
%
%   P = propagateParticleMultiMaterial(matArray,geometry,P,T) generalizes
%   propagateParticleSmallDt.m to a domain built from several homogeneous
%   regions (matArray{1}, matArray{2}, ...), with the interfaces between
%   them given genuine Zoeppritz transmission instead of the single-medium
%   code's solid-air reflection. A particle carries which region it is
%   currently in as P.matIdx (1-based, into matArray); each sub-step all
%   particles advance using THEIR OWN region's velocity and mean free
%   time, and P/SV/SH are reflected or transmitted (with mode conversion
%   on either branch) according to MaterialInterfaceTwoSolid or (when one
%   side has vs==0) MaterialInterfaceSolidFluid.
%
%   Interface boundaries are geometry.bnd entries with a .mat = [i j]
%   field naming the two region indices they separate, and dir == 1, 2 or
%   3 (a plane normal to x, y or z respectively: layered media, matching
%   the axis convention propagateParticleSmallDt.m already uses for its
%   own dir<4 boundaries), dir == 4 (a cylinder: concentric regions, x-y
%   radius) or dir == 5 (a sphere: concentric shells, full 3D radius --
%   an Earth-like core/mantle/crust model). A particle whose matIdx is i
%   or j is tested against it, all others ignore it. Boundary entries
%   without .mat behave exactly as in the single-material code
%   (reflective/absorbing against vacuum, using whichever material the
%   particle currently carries).
%
%   dir == 5 is specific to THIS function: propagateParticleSmallDt.m has
%   no spherical boundary at all (only dir<4 planes and dir==4 cylindrical
%   radius), but this function's "ordinary boundary" handling never calls
%   into it -- it uses its own ordinaryBoundaryDecision (below), a
%   per-particle-normal generalisation of MaterialInterface.m's P/SV/SH
%   split kept local to this file, so no upstream change was needed to add
%   this geometry.
%
%   This is intentionally a separate function rather than a rewrite of
%   propagateParticleSmallDt.m: the validated single-material path stays
%   untouched.
%
%   See also MATERIALINTERFACETWOSOLID, MATERIALINTERFACESOLIDFLUID,
%   PROPAGATEPARTICLESMALLDT.

T = T - mean(P.t);
if T <= 0
    return
end

Nmat = numel(matArray);
isFluid = cellfun(@(m) m.vs==0, matArray);

% Sub-step, small against every material's own mean free time AND against
% the domain geometry: a step must not be able to cross more than one
% boundary at once, since only a single crossing is detected and resolved
% per boundary per sub-step -- true along each cartesian axis separately,
% and equally true radially (the gap between two concentric cylindrical
% interfaces must not be crossable in one step either).
mft = cellfun(@(m) min(m.meanFreeTime(:)), matArray);
dt = 0.01*min(mft);

vmax = max(cellfun(@(m) max(m.vp,m.vs), matArray));
gaps = [];
for axis_ = 1:3
    avals = arrayfun(@(b) b.val, geometry.bnd([geometry.bnd.dir]==axis_));
    if numel(avals) > 1, gaps(end+1) = min(diff(sort(avals))); end %#ok<AGROW>
end
rvals = arrayfun(@(b) b.val, geometry.bnd([geometry.bnd.dir]==4));
if numel(rvals) > 1, gaps(end+1) = min(diff(sort(rvals))); end %#ok<AGROW>
Rvals = arrayfun(@(b) b.val, geometry.bnd([geometry.bnd.dir]==5));
if numel(Rvals) > 1, gaps(end+1) = min(diff(sort(Rvals))); end %#ok<AGROW>
% For a radial (dir 4 or 5) boundary, the distance from the axis/origin
% along a straight chord is not monotonic -- it dips to a minimum at the
% point of closest approach. A substep long enough to cross deep into a
% SMALL shell and back out within itself has both endpoints on the same
% side of that shell's radius, so the endpoint sign-change test below
% never sees the crossing at all (found by this session's own spherical
% validation: a small but nonzero leak at the smallest, innermost shell).
% Capping the step length against the smallest radius present too (not
% just the gaps between radii) keeps that miss-probability negligible.
allR = [rvals, Rvals];
allR = allR(allR > 0);
if ~isempty(allR), gaps(end+1) = min(allR); end %#ok<AGROW>
if ~isempty(gaps)
    gap = min(gaps);
    % 0.0971, not a round 0.1: an exact fraction here combined with round
    % problem geometry (radii like 0.10/0.25/0.40) is exactly what causes
    % the near-degenerate case the tooClose guard above exists for --
    % this just makes that coincidence far less likely to begin with.
    if gap > 0 && isfinite(gap), dt = min(dt, 0.0971*gap/vmax); end
end

dt = min(dt, T);
Nt = ceil(T/dt);
dt = T/Nt;

% Zoeppritz table caches, one per unordered material pair actually present
% at an interface: either a two-solid table, or -- when either side has
% vs==0 -- a solid-incident pair of tables plus a fluid-incident one (see
% MaterialInterfaceSolidFluid.m for why a solid-fluid pair needs both).
twoSolidCache = containers.Map('KeyType','char','ValueType','any');
solidFluidCache = containers.Map('KeyType','char','ValueType','any');
ordinaryCache = cell(1,Nmat); ordinaryCache2 = cell(1,Nmat); % per-material Zoeppritz-vs-air interpolants (Rpp, Rsvsv), built lazily
j1_deg = linspace(0,90,361);
pairKey = @(i,j) sprintf('%d_%d', min(i,j), max(i,j));

for i2 = 1:numel(geometry.bnd)
    bnd = geometry.bnd(i2);
    if bnd.dir < 1 || bnd.dir > 5
        error('propagateParticleMultiMaterial only supports boundary dir 1-3 (a plane normal to x, y or z), 4 (a cylinder) or 5 (a sphere).');
    end
    if ~(isfield(bnd,'mat') && ~isempty(bnd.mat)), continue; end
    i = bnd.mat(1); j = bnd.mat(2);
    key = pairKey(i,j);
    if isFluid(i) || isFluid(j)
        if isFluid(i) && isFluid(j), continue; end % fluid-fluid: closed-form R(theta), no table (fluidFluidDecision)
        if isKey(solidFluidCache, key), continue; end
        if isFluid(j), solidIdx = i; fluidIdx = j; else, solidIdx = j; fluidIdx = i; end
        mS = matArray{solidIdx}; mF = matArray{fluidIdx};
        zS.P  = MaterialClass.ZoeppritzFluid(j1_deg, mS.vp,mS.vs,mS.rho, mF.vp,0,mF.rho);
        zS.SV = zS.P; % the same call already carries both incidence blocks
        zS.j1_deg = j1_deg; zS.vpFluid = mF.vp;
        zF = MaterialClass.ZoeppritzFluidIncident(j1_deg, mS.vp,mS.vs,mS.rho, mF.vp,mF.rho);
        solidFluidCache(key) = struct('solidIdx',solidIdx,'fluidIdx',fluidIdx,'zSolid',zS,'zFluid',zF); %#ok<NASGU>
    else
        if isKey(twoSolidCache,key), continue; end
        mA = matArray{i}; mB = matArray{j};
        out = MaterialClass.ZoeppritzSolid(j1_deg,mA.vp,mA.vs,mA.rho,mB.vp,mB.vs,mB.rho);
        out.j1_deg = j1_deg;
        twoSolidCache(key) = out; %#ok<NASGU>
    end
end

if ~isfield(P,'matIdx')
    P.matIdx = ones(size(P.p));
end

for i1 = 1:Nt

    alive = P.alive;
    if ~any(alive), break; end

    x1 = P.x;

    % straight-line advance, each particle at its own region's velocity
    for k = 1:Nmat
        m = alive & P.matIdx==k;
        if ~any(m), continue; end
        mat = matArray{k};
        p = m & P.p; s = m & ~P.p;
        P.x(p,:) = P.x(p,:) + (mat.vp*dt).*P.dir(p,:);
        if any(s), P.x(s,:) = P.x(s,:) + (mat.vs*dt).*P.dir(s,:); end
    end

    % --- boundaries ---
    for i2 = 1:numel(geometry.bnd)
        bnd = geometry.bnd(i2);

        if bnd.dir < 4
            r0 = x1(:,bnd.dir); r1v = P.x(:,bnd.dir);
        elseif bnd.dir == 4 % radial (cylinder, x-y plane)
            r0 = hypot(x1(:,1),x1(:,2));
            r1v = hypot(P.x(:,1),P.x(:,2));
        else % dir == 5, radial (sphere, full 3D)
            r0 = vecnorm(x1,2,2);
            r1v = vecnorm(P.x,2,2);
        end
        % A radial (dir 4/5) step length that happens to be an exact
        % multiple of a shell's radius -- an everyday coincidence with
        % round problem geometry (R = 0.10/0.25/0.40 and the geometric dt
        % cap above combine to exactly that) -- lands a purely radial
        % particle EXACTLY on the shell. The very next substep then has
        % r0-bnd.val sitting at plain floating-point noise around zero,
        % and whether that noise is a few ULPs positive or negative --
        % effectively arbitrary -- decides whether the endpoint sign-change
        % test below sees a crossing at all.
        %
        % An EARLIER version of this guard SKIPPED the test for a tooClose
        % particle, meaning to let it "resolve itself" next substep. That
        % is wrong and was found to corrupt trajectories (a particle
        % genuinely transmitting outward from r0==bnd.val to r1v>bnd.val
        % has BOTH endpoints reported on the outside from then on, so the
        % sign-change test can never see the crossing again -- the event
        % is not deferred, it is lost, and the particle sails straight
        % through every later boundary undetected). Found by this
        % session's own 3-shell spherical validation (MCRTE/todo.rtf).
        %
        % The fix is the opposite: ALWAYS resolve a tooClose particle as
        % crossed, right now. This is safe because the boundary decision
        % itself (MaterialInterfaceTwoSolid / MaterialInterfaceSolidFluid /
        % ordinaryBoundaryDecision) determines reflect vs. transmit from
        % the particle's ACTUAL DIRECTION at the interface, not from which
        % side its position nominally reports -- a particle sitting right
        % at the boundary is exactly a normal-incidence hit, and the
        % Zoeppritz decision at theta=0 is well defined regardless of
        % which side floating-point rounding puts r0 on.
        tooClose = abs(r0-bnd.val) < 1e-9*max(1,abs(bnd.val));
        crossed = alive & (tooClose | ((r1v-bnd.val).*(r0-bnd.val) < 0));
        if ~any(crossed), continue; end

        nrm = localNormal(P, bnd);            % P.N x 3, only rows in `crossed` are used

        isInterface = isfield(bnd,'mat') && ~isempty(bnd.mat);

        if ~isInterface
            bndType = 'reflective';
            if isfield(bnd,'type') && ~isempty(bnd.type), bndType = lower(bnd.type); end
            if strcmp(bndType,'absorbing')
                P.alive(crossed) = false;
                continue
            end
            P = foldPosition(P, crossed, bnd);
            for k = 1:Nmat
                ind = crossed & P.matIdx==k;
                if ~any(ind), continue; end
                if isFluid(k)
                    % A fluid against vacuum has no shear modulus on
                    % either side of the (non-)interface: it is not a
                    % milder version of the solid-vs-air case, it is
                    % degenerate for it (getZoeppritzCached ultimately
                    % calls ZoeppritzFluid(fluid,...,air,...), two
                    % zero-shear media at once, which produced a singular
                    % boundary-condition matrix when this was tried).
                    % Physically a fluid free surface simply reflects the
                    % P wave specularly (there is nothing else it CAN do:
                    % no S to convert into, no shear traction to speak
                    % of), so that is applied directly, no Zoeppritz call.
                    n_idx = nrm(ind,:);
                    dn = sum(P.dir(ind,:).*n_idx,2);
                    P.dir(ind,:) = P.dir(ind,:) - 2*dn.*n_idx;
                    % same post-reflection nudge as setDirection, and for
                    % the same reason: fold(val)=val is a no-op sitting
                    % exactly on the boundary, so this stops the very next
                    % substep's tooClose check reprocessing this event.
                    P.x(ind,:) = P.x(ind,:) + (1e-9*max(1,vecnorm(P.x(ind,:),2,2))) .* P.dir(ind,:);
                    continue
                end
                if isempty(ordinaryCache{k})
                    [outZ,~] = matArray{k}.getZoeppritzCached();
                    ordinaryCache{k} = griddedInterpolant(outZ.j1_deg, outZ.E_Rpp, 'linear', 'nearest');
                    ordinaryCache2{k} = griddedInterpolant(outZ.j1_deg, outZ.E_Rsvsv, 'linear', 'nearest'); %#ok<AGROW>
                end
                o = ordinaryBoundaryDecision(matArray{k}, ordinaryCache{k}, ordinaryCache2{k}, P, nrm, ind);
                P = applyInterfaceDecision(P, o, nrm, matArray{k}.vp, matArray{k}.vs, ...
                    matArray{k}.vp, matArray{k}.vs, k, k, []);
            end
            continue
        end

        i = bnd.mat(1); j = bnd.mat(2);
        indI = crossed & P.matIdx==i;
        indJ = crossed & P.matIdx==j;
        if ~any(indI) && ~any(indJ), continue; end

        key = pairKey(i,j);
        if isFluid(i) && isFluid(j)
            matI = matArray{i}; matJ = matArray{j};
            if any(indI)
                o = fluidFluidDecision(matI, matJ, P, nrm, indI);
                P = applyInterfaceDecision(P, o, nrm, matI.vp, 0, matJ.vp, 0, i, j, bnd);
            end
            if any(indJ)
                o = fluidFluidDecision(matJ, matI, P, -nrm, indJ);
                P = applyInterfaceDecision(P, o, -nrm, matJ.vp, 0, matI.vp, 0, j, i, bnd);
            end
        elseif isFluid(i) || isFluid(j)
            sf = solidFluidCache(key);
            matI = matArray{i}; matJ = matArray{j};
            if any(indI)
                fromFluid = (i == sf.fluidIdx);
                tbl = ternary(fromFluid, sf.zFluid, sf.zSolid);
                o = matI.MaterialInterfaceSolidFluid(fromFluid, tbl, P, nrm, indI);
                vsI = ternary(fromFluid, 0, matI.vs);
                P = applyInterfaceDecision(P, o, nrm, matI.vp, vsI, matJ.vp, matJ.vs, i, j, bnd);
            end
            if any(indJ)
                fromFluid = (j == sf.fluidIdx);
                tbl = ternary(fromFluid, sf.zFluid, sf.zSolid);
                o = matJ.MaterialInterfaceSolidFluid(fromFluid, tbl, P, -nrm, indJ);
                vsJ = ternary(fromFluid, 0, matJ.vs);
                P = applyInterfaceDecision(P, o, -nrm, matJ.vp, vsJ, matI.vp, matI.vs, j, i, bnd);
            end
        else
            out = twoSolidCache(key);
            matI = matArray{i}; matJ = matArray{j};
            if any(indI)
                o = matI.MaterialInterfaceTwoSolid(matJ, out, P, nrm, indI);
                P = applyInterfaceDecision(P, o, nrm, matI.vp, matI.vs, matJ.vp, matJ.vs, i, j, bnd);
            end
            if any(indJ)
                o = matJ.MaterialInterfaceTwoSolid(matI, out, P, -nrm, indJ);
                P = applyInterfaceDecision(P, o, -nrm, matJ.vp, matJ.vs, matI.vp, matI.vs, j, i, bnd);
            end
        end
    end

    if ~any(P.alive), break; end

    % --- scattering, per material ---
    for k = 1:Nmat
        alive = P.alive;
        m = alive & P.matIdx==k;
        if ~any(m), continue; end
        mat = matArray{k};
        p = m & P.p; s = m & ~P.p;

        lambda = zeros(P.N,1);
        lambda(p) = dt/mat.meanFreeTime(1);
        if ~isFluid(k), lambda(s) = dt/mat.meanFreeTime(2); end
        scat = rand(P.N,1) < lambda; % single-event approx, dt << meanFreeTime
        scat = scat & m;
        Nscat = nnz(scat);
        if Nscat == 0, continue; end

        if P.d == 2 % in-plane: perp only flips side (as propagateParticleSmallDt)
            flip = rand(Nscat,1) < 0.5;
            P.perp(scat,:) = (1 - 2*flip).*P.perp(scat,:);
        else
            phi = 2*pi*rand(Nscat,1);
            P.perp(scat,:) = cos(phi).*P.perp(scat,:) + sin(phi).*cross(P.dir(scat,:),P.perp(scat,:));
        end

        pS = scat & p; sS = scat & s;
        if ~isFluid(k)
            P.p( pS & (rand(P.N,1)>mat.P2P) ) = false;
            P.p( sS & (rand(P.N,1)>mat.S2S) ) = true;
        end % a fluid never converts mode: P.p stays true throughout

        theta = zeros(P.N,1);
        u = rand(P.N,1);
        if iscell(mat.invcdf)
            theta(  P.p & pS ) = mat.invcdf{1,1}( u(  P.p & pS ) );
        else % acoustic MaterialClass: a single invcdf handle
            theta(  P.p & pS ) = mat.invcdf( u(  P.p & pS ) );
        end
        if ~isFluid(k)
            theta( ~P.p & pS ) = mat.invcdf{1,2}( u( ~P.p & pS ) );
            theta(  P.p & sS ) = mat.invcdf{2,1}( u(  P.p & sS ) );
            theta( ~P.p & sS ) = mat.invcdf{2,2}( u( ~P.p & sS ) );
        end

        dir = P.dir;
        P.dir(scat,:)  =  cos(theta(scat)).*dir(scat,:) + sin(theta(scat)).*P.perp(scat,:);
        P.dir(scat,:)  = P.dir(scat,:)./sqrt(sum(P.dir(scat,:).^2,2));
        P.perp(scat,:) = -sin(theta(scat)).*dir(scat,:) + cos(theta(scat)).*P.perp(scat,:);
        P.perp(scat,:) = P.perp(scat,:)./sqrt(sum(P.perp(scat,:).^2,2));
    end
end

P.t = P.t + T;
end

% -------------------------------------------------------------------------
function v = ternary(cond, a, b)
if cond, v = a; else, v = b; end
end

% -------------------------------------------------------------------------
function o = fluidFluidDecision(matA, matB, P, nrm, ind)
% Acoustic A -> B interface: Snell sin2 = cB/cA sin1, energy reflection
% R = ((ZB cos1 - ZA cos2)/(ZB cos1 + ZA cos2))^2, Z = rho*c, total
% reflection beyond the critical angle; reflect with probability R.
cos1 = min(1, abs(sum(P.dir(ind,:).*nrm(ind,:), 2)));
sin2 = (matB.vp/matA.vp) * sqrt(1 - cos1.^2);
R = ones(size(cos1));
ok = sin2 < 1;
cos2 = sqrt(1 - sin2(ok).^2);
if isempty(matA.rho) || isempty(matB.rho)
    error('Fluid-fluid interface: set rho on both materials (impedance Z = rho*c).');
end
ZA = matA.rho*matA.vp; ZB = matB.rho*matB.vp;
R(ok) = ((ZB*cos1(ok) - ZA*cos2) ./ (ZB*cos1(ok) + ZA*cos2)).^2;
refl = rand(size(R)) <= R;
o.pparticle = ind;
o.p2p  = refl;
o.p2Tp = ~refl;
end

% -------------------------------------------------------------------------
function nrm = localNormal(P, bnd)
% outward normal at every particle's CURRENT position, pointing from the
% side named bnd.mat(1) [or, for an ordinary boundary, from inside the
% domain] towards bnd.mat(2)/outside. Only meaningful (and only used) at
% the rows that actually crossed this boundary this sub-step.
nrm = zeros(P.N,3);
if bnd.dir < 4
    nrm(:,bnd.dir) = 1;
elseif bnd.dir == 4 % radial (cylinder, x-y plane)
    r = hypot(P.x(:,1),P.x(:,2));
    r = max(r,eps);
    nrm(:,1) = P.x(:,1)./r;
    nrm(:,2) = P.x(:,2)./r;
else % dir == 5, radial (sphere, full 3D)
    r = vecnorm(P.x,2,2);
    r = max(r,eps);
    nrm = P.x ./ r;
end
end

% -------------------------------------------------------------------------
function P = foldPosition(P, idx, bnd)
% mirror a particle's position back across the boundary it just crossed
% (used for reflection, and for the near-critical transmission fallback);
% never called for a genuinely transmitted particle, which keeps its
% already-overshot position.
if ~any(idx), return; end
if bnd.dir < 4
    P.x(idx,bnd.dir) = 2*bnd.val - P.x(idx,bnd.dir);
elseif bnd.dir == 4
    r = hypot(P.x(idx,1),P.x(idx,2));
    newR = 2*bnd.val - r;
    scale = newR ./ max(r,eps);
    P.x(idx,1) = P.x(idx,1).*scale;
    P.x(idx,2) = P.x(idx,2).*scale;
else % dir == 5, radial (sphere, full 3D): same fold, all three components
    r = vecnorm(P.x(idx,:),2,2);
    newR = 2*bnd.val - r;
    scale = newR ./ max(r,eps);
    P.x(idx,:) = P.x(idx,:).*scale;
end
end

% -------------------------------------------------------------------------
function o = ordinaryBoundaryDecision(mat, RppInterp, RsvsvInterp, P, nrm, ind)
% Mirrors MaterialInterface.m's P/SV/SH split and Monte Carlo reflection
% decision (solid against vacuum/air), generalised to a per-particle
% normal so it also works at a cylindrical (dir==4) ordinary boundary,
% where MaterialInterface.m's own repmat(n,nP,1) (a single shared normal)
% does not apply. Same physics, same cached Zoeppritz-vs-air tables
% (mat.getZoeppritzCached), only the normal's source differs.
pIdx =  P.p  & ind;
sIdx = ~P.p  & ind;

partP_dir = P.dir(pIdx,:);
nP = size(partP_dir,1);
if nP > 0
    nRepP = nrm(pIdx,:);
    incAngP = acosd(min(1,abs(dot(partP_dir, nRepP, 2))));
    Rpp_coeff = RppInterp(incAngP);
    p2p  = rand(nP,1) <= Rpp_coeff;
    p2sv = ~p2p;
else
    p2p = []; p2sv = [];
end

partS_dir  = P.dir(sIdx,:);
partS_perp = P.perp(sIdx,:);
nS = size(partS_dir,1);
if nS > 0
    nRepS = nrm(sIdx,:);
    eSH = cross(partS_dir, nRepS, 2);
    norm_eSH = vecnorm(eSH, 2, 2);
    nearNormal = norm_eSH < eps;
    if any(nearNormal)
        nNear = nRepS(nearNormal,:);
        [~,referenceIndex] = min(abs(nNear),[],2);
        t = zeros(size(nNear));
        t(sub2ind(size(t),(1:size(t,1))',referenceIndex)) = 1;
        t = t - sum(t.*nNear,2).*nNear;
        eSH(nearNormal,:)    = t;
        norm_eSH(nearNormal) = vecnorm(t, 2, 2);
    end
    eSH = eSH ./ norm_eSH;
    eSV = cross(eSH, partS_dir, 2);
    eSV = eSV ./ vecnorm(eSV, 2, 2);

    aSV = dot(partS_perp, eSV, 2);
    aSH = dot(partS_perp, eSH, 2);
    E_tot = aSV.^2 + aSH.^2 + eps;
    isSV = rand(size(aSV)) < (aSV.^2 ./ E_tot);

    incAngS = acosd(min(1,abs(dot(partS_dir, nRepS, 2))));
    incAngSV = incAngS(isSV);

    nSV = nnz(isSV);
    if nSV > 0
        Rsvsv_coeff = RsvsvInterp(incAngSV);
        sv2sv = rand(nSV,1) <= Rsvsv_coeff;
        sv2p  = ~sv2sv;
    else
        sv2sv = []; sv2p = [];
    end
    sh2sh = true(nnz(~isSV),1);
else
    isSV = false(0,1);
    sv2sv = []; sv2p = []; sh2sh = [];
end

o.pparticle  = pIdx;
o.sparticle  = sIdx;
o.svparticle = isSV;
o.shparticle = ~isSV;
o.p2p = p2p; o.p2sv = p2sv;
o.sv2sv = sv2sv; o.sv2p = sv2p;
o.sh2sh = sh2sh;
end

% -------------------------------------------------------------------------
function P = applyInterfaceDecision(P, o, nrm, vpA, vsA, vpB, vsB, idxA, idxB, bnd)
% Apply a MaterialInterface*-style decision struct to P.x/P.dir/P.perp/
% P.p/P.matIdx. nrm is the per-particle outward normal (P.N x 3, pointing
% from medium A into medium B); reflection folds the particle's position
% back across bnd and keeps it in medium A, transmission (fields ending
% in T) moves it into medium B with Snell's law, keeping its already-
% overshot position. bnd is [] for the ordinary (never-transmits) case,
% where position was already folded in bulk by the caller.
%
% Every branch is applied via a direct setDirection(P, ...) call, not an
% anonymous-function wrapper: an anonymous function captures P BY VALUE at
% the point it is defined, so creating it once and reassigning P on each
% call would make every subsequent call see that first, stale P and undo
% the previous branches' updates.

if isfield(o,'p2p') && any(o.p2p)
    idx = o.pparticle; idx(idx) = o.p2p;
    P = setDirection(P, idx, nrm, vpA, vpA, false, 'same', idxA, bnd);
end
if isfield(o,'p2sv') && any(o.p2sv)
    idx = o.pparticle; idx(idx) = o.p2sv;
    P = setDirection(P, idx, nrm, vpA, vsA, false, 'conv', idxA, bnd);
end
if isfield(o,'p2Tp') && any(o.p2Tp)
    idx = o.pparticle; idx(idx) = o.p2Tp;
    P = setDirection(P, idx, nrm, vpA, vpB, true, 'same', idxB, bnd);
end
if isfield(o,'p2Tsv') && any(o.p2Tsv)
    idx = o.pparticle; idx(idx) = o.p2Tsv;
    P = setDirection(P, idx, nrm, vpA, vsB, true, 'conv', idxB, bnd);
end

if isfield(o,'sv2sv') && any(o.sv2sv)
    isSV = false(P.N,1); isSV(o.sparticle) = o.svparticle;
    idx = false(P.N,1); idx(isSV) = o.sv2sv;
    P = setDirection(P, idx, nrm, vsA, vsA, false, 'same', idxA, bnd);
end
if isfield(o,'sv2p') && any(o.sv2p)
    isSV = false(P.N,1); isSV(o.sparticle) = o.svparticle;
    idx = false(P.N,1); idx(isSV) = o.sv2p;
    P = setDirection(P, idx, nrm, vsA, vpA, false, 'conv', idxA, bnd);
end
if isfield(o,'sv2Tsv') && any(o.sv2Tsv)
    isSV = false(P.N,1); isSV(o.sparticle) = o.svparticle;
    idx = false(P.N,1); idx(isSV) = o.sv2Tsv;
    P = setDirection(P, idx, nrm, vsA, vsB, true, 'same', idxB, bnd);
end
if isfield(o,'sv2Tp') && any(o.sv2Tp)
    isSV = false(P.N,1); isSV(o.sparticle) = o.svparticle;
    idx = false(P.N,1); idx(isSV) = o.sv2Tp;
    P = setDirection(P, idx, nrm, vsA, vpB, true, 'conv', idxB, bnd);
end

if isfield(o,'sh2sh') && any(o.sh2sh)
    isSH = false(P.N,1); isSH(o.sparticle) = o.shparticle;
    idx = false(P.N,1); idx(isSH) = o.sh2sh;
    P = setDirection(P, idx, nrm, vsA, vsA, false, 'same', idxA, bnd);
end
if isfield(o,'sh2Tsh') && any(o.sh2Tsh)
    isSH = false(P.N,1); isSH(o.sparticle) = o.shparticle;
    idx = false(P.N,1); idx(isSH) = o.sh2Tsh;
    P = setDirection(P, idx, nrm, vsA, vsB, true, 'shT', idxB, bnd);
end
end

% -------------------------------------------------------------------------
function P = setDirection(P, idx, nrm, vin, vout, transmit, mode, newMatIdx, bnd)
if ~any(idx), return; end

n_idx = nrm(idx,:);
d_old = P.dir(idx,:);
d_n = sum(d_old.*n_idx,2);
d_t = d_old - d_n.*n_idx;

ratio = vout/vin;
d_new_t = ratio*d_t;
sin2 = sum(d_new_t.^2,2);

valid = sin2 <= 1;
if transmit && any(~valid)
    % ZoeppritzSolid/Fluid already drove E_T to ~0 there, so this branch
    % should be rare/empty; if hit anyway (angle-grid interpolation slop
    % right at the critical angle), fall back to specular reflection
    % rather than producing a NaN direction.
    fb = idx; fb(idx) = ~valid;
    if ~isempty(bnd), P = foldPosition(P, fb, bnd); end
    dn_fb = sum(P.dir(fb,:).*nrm(fb,:),2);
    P.dir(fb,:) = P.dir(fb,:) - 2*dn_fb.*nrm(fb,:);
    dp_fb = sum(P.perp(fb,:).*P.dir(fb,:),2);
    P.perp(fb,:) = P.perp(fb,:) - dp_fb.*P.dir(fb,:);
    P.perp(fb,:) = P.perp(fb,:) ./ vecnorm(P.perp(fb,:),2,2);
    P.x(fb,:) = P.x(fb,:) + (1e-9*max(1,vecnorm(P.x(fb,:),2,2))) .* P.dir(fb,:); % same post-reflection nudge
    sub = find(idx); idx = false(size(idx)); idx(sub(valid)) = true;
    n_idx = n_idx(valid,:); d_old = d_old(valid,:); d_n = d_n(valid); %#ok<NASGU>
    d_new_t = d_new_t(valid,:); sin2 = sin2(valid);
    if ~any(idx), return; end
end

if ~transmit && ~isempty(bnd)
    P = foldPosition(P, idx, bnd);
end

cos_val = sqrt(max(0,1-sin2));
sgn = -sign(d_n); if transmit, sgn = sign(d_n); end
d_new = d_new_t + sgn.*cos_val.*n_idx;
d_new = d_new ./ vecnorm(d_new,2,2);
P.dir(idx,:) = d_new;

% Nudge the position a hair along the freshly resolved direction. A
% particle that was sitting exactly at the boundary before this event
% (r0 == bnd.val, the tooClose case above) is STILL sitting exactly there
% after a reflection: folding r == val is a no-op. Without this nudge, the
% very next substep sees the same r0 == bnd.val again and re-triggers the
% tooClose branch -- reprocessing the SAME event a second time, now with
% the just-reflected (opposite) direction, which can talk it into a wrong
% second decision (observed this session: a particle reflects correctly,
% then gets "transmitted" on the immediately following substep purely
% because it never numerically left the boundary). The nudge is along the
% direction actually just computed, not a guess: it is now the one piece
% of unambiguous information available about which way this particle is
% really going.
P.x(idx,:) = P.x(idx,:) + (1e-9*max(1,vecnorm(P.x(idx,:),2,2))) .* d_new;

switch mode
    case 'same'
        % Gram-Schmidt against the new direction rather than a sign flip:
        % a sign flip only keeps perp orthogonal for reflection (where
        % dir's normal component flips too), not for transmission (where
        % it does not).
        P.perp(idx,:) = P.perp(idx,:) - sum(P.perp(idx,:).*P.dir(idx,:),2).*P.dir(idx,:);
        P.perp(idx,:) = P.perp(idx,:) ./ vecnorm(P.perp(idx,:),2,2);
    case {'conv','shT'}
        n_use = nrm(idx,:);
        axis_rot = cross(P.dir(idx,:), n_use, 2);
        na = vecnorm(axis_rot,2,2);
        bad = na < 1e-9;
        if any(bad)
            alt = repmat([1 0 0], nnz(bad), 1);
            sub = find(idx); axis_rot(bad,:) = cross(P.dir(sub(bad),:), alt, 2);
            na(bad) = vecnorm(axis_rot(bad,:),2,2);
        end
        axis_rot = axis_rot ./ na;
        if strcmp(mode,'shT')
            % A transmitted SH wave stays polarized OUT of the plane of
            % incidence, and axis_rot = cross(dir,n) is that plane's normal
            % -- already a unit vector orthogonal to the new direction, so it
            % is the polarization itself.  cross(axis_rot,dir) is the SV
            % basis vector: right for 'conv' (P<->SV really does come out
            % SV-polarized), wrong here, where it silently turns every
            % transmitted SH wave into an SV one.  The channel probability
            % was always correct; only the outgoing polarization was not,
            % which is why the direct MaterialInterface* tests could not see
            % it -- they read the routine's own channel flags and never look
            % at P.perp afterwards.
            P.perp(idx,:) = axis_rot;
        else
            P.perp(idx,:) = cross(axis_rot, P.dir(idx,:), 2);
        end
end

if strcmp(mode,'conv')
    P.p(idx) = ~P.p(idx);
end
P.matIdx(idx) = newMatIdx;
end
