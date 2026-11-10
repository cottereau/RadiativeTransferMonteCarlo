function o = MaterialInterfaceTwoSolid(matA, matB, zoepAB, P, n, ind)
%MATERIALINTERFACETWOSOLID  Reflect/transmit particles at a solid-solid weld.
%
%   o = MaterialInterfaceTwoSolid(matA, matB, zoepAB, P, n, ind)
%   mirrors MaterialInterface (which assumes a solid-air free surface, pure
%   reflection) but for a welded interface between two genuine solids: for
%   each incident P/SV/SH particle it Monte-Carlo selects one of four
%   outcomes -- reflected same mode, reflected converted mode, transmitted
%   same mode, transmitted converted mode -- weighted by ZoeppritzSolid's
%   energy coefficients (E_Rpp,E_Rpsv,E_Tpp,E_Tpsv for P incidence;
%   E_Rsvsv,E_Rsp,E_Tsvsv,E_Tsp for SV; E_Rsh,E_Tsh for SH). Reflection past
%   the critical angle needs no special case: ZoeppritzSolid's complex-angle
%   formulation already drives E_T to exactly zero there (real(cos(evanescent
%   angle)) = 0), so the remaining probability mass is automatically all
%   reflection.
%
%   matA is the medium the incident particles are currently in; matB is the
%   medium across the interface. zoepAB is the table from
%   MaterialClass.ZoeppritzSolid(j1_deg, matA.vp,matA.vs,matA.rho,
%   matB.vp,matB.vs,matB.rho), passed in (rather than recomputed here) so
%   the caller can build it once per material pair instead of once per
%   sub-step -- it depends only on the four velocities and two densities,
%   not on any particle state.
%
%   P   - particle structure (fields: p, dir, perp, ...)
%   n   - interface normal (1x3), pointing from medium A into medium B
%   ind - logical index of particles hitting this interface (already
%         confirmed to be in medium A: P.matIdx(ind) == the index of matA)
%
%   Output o carries, for each of P/SV/SH incidence, four logical masks
%   (same/converted x reflect/transmit) plus the reflection/transmission
%   angle (degrees, from the interface normal) for each -- everything the
%   caller needs to update P.dir/P.perp/P.matIdx, exactly as
%   propagateParticleSmallDt.m already does with MaterialInterface's output
%   for the free-surface case. This function makes the decision; it does
%   not touch P itself.
%
%   See also MATERIALINTERFACE, ZOEPPRITZSOLID, PROPAGATEPARTICLEMULTIMATERIAL.

if isempty(matA.rho) || isempty(matB.rho)
    error(['Please define the density of both materials before calling ' ...
        'MaterialInterfaceTwoSolid.']);
end

j1_deg = zoepAB.j1_deg;
lookup = @(data, ang) interp1(j1_deg, data, ang, 'linear', 'extrap');

vpA = matA.vp; vsA = matA.vs;
vpB = matB.vp; vsB = matB.vs;

% Snell's law angles, mirroring the angles.* handles in MaterialClass.ZoeppritzVsAir
angRpp  = @(th) th;                        % specular
angRpsv = @(th) asind(min(1,(sind(th)/vpA)*vsA));
angTpp  = @(th) asind(min(1,(sind(th)/vpA)*vpB));
angTpsv = @(th) asind(min(1,(sind(th)/vpA)*vsB));

angRsvsv = @(th) th;
angRsvp  = @(th) asind(min(1,(sind(th)/vsA)*vpA));
angTsvsv = @(th) asind(min(1,(sind(th)/vsA)*vsB));
angTsvp  = @(th) asind(min(1,(sind(th)/vsA)*vpB));

angTsh = @(th) asind(min(1,(sind(th)/vsA)*vsB));

% --- split P vs S incident particles ---
pIdx =  P.p  & ind;
sIdx = ~P.p  & ind;

partP_dir = P.dir(pIdx,:);
nP = size(partP_dir,1);
if nP > 0
    nRepP = n(pIdx,:);
    incAngP = acosd(min(1,abs(dot(partP_dir, nRepP, 2))));
else
    incAngP = [];
end

partS_dir  = P.dir(sIdx,:);
partS_perp = P.perp(sIdx,:);
nS = size(partS_dir,1);
if nS > 0
    nRepS = n(sIdx,:);
    eSH = cross(partS_dir, nRepS, 2);
    norm_eSH = vecnorm(eSH, 2, 2);
    nearNormal = norm_eSH < eps;
    if any(nearNormal)
        % n is now a per-particle normal array (needed for the cylindrical
        % boundary, where it genuinely varies row to row): build a
        % reference axis per near-normal row instead of one shared axis.
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
    incAngSH = incAngS(~isSV);
else
    isSV = false(0,1);
    incAngSV = []; incAngSH = [];
end

% --- P incidence: reflect-P / reflect-SV / transmit-P / transmit-SV ---
if nP > 0
    Rpp  = lookup(zoepAB.E_Rpp,  incAngP);
    Rpsv = lookup(zoepAB.E_Rpsv, incAngP);
    Tpp  = lookup(zoepAB.E_Tpp,  incAngP);
    Tpsv = lookup(zoepAB.E_Tpsv, incAngP);
    Wtot = Rpp + Rpsv + Tpp + Tpsv;          % ~1, renormalise against roundoff
    u = rand(nP,1) .* Wtot;
    c1 = Rpp; c2 = c1+Rpsv; c3 = c2+Tpp;
    p2p   = u <= c1;
    p2sv  = u > c1  & u <= c2;
    p2Tp  = u > c2  & u <= c3;
    p2Tsv = u > c3;

    angP2P   = angRpp(incAngP(p2p));
    angP2SV  = angRpsv(incAngP(p2sv));
    angP2Tp  = angTpp(incAngP(p2Tp));
    angP2Tsv = angTpsv(incAngP(p2Tsv));
else
    p2p=[]; p2sv=[]; p2Tp=[]; p2Tsv=[];
    angP2P=[]; angP2SV=[]; angP2Tp=[]; angP2Tsv=[];
end

% --- SV incidence: reflect-SV / reflect-P / transmit-SV / transmit-P ---
nSV = numel(incAngSV);
if nSV > 0
    Rsvsv = lookup(zoepAB.E_Rsvsv, incAngSV);
    Rsp   = lookup(zoepAB.E_Rsp,   incAngSV);
    Tsvsv = lookup(zoepAB.E_Tsvsv, incAngSV);
    Tsp   = lookup(zoepAB.E_Tsp,   incAngSV);
    Wtot = Rsvsv + Rsp + Tsvsv + Tsp;
    u = rand(nSV,1) .* Wtot;
    c1 = Rsvsv; c2 = c1+Rsp; c3 = c2+Tsvsv;
    sv2sv  = u <= c1;
    sv2p   = u > c1 & u <= c2;
    sv2Tsv = u > c2 & u <= c3;
    sv2Tp  = u > c3;

    angSV2SV  = angRsvsv(incAngSV(sv2sv));
    angSV2P   = angRsvp(incAngSV(sv2p));
    angSV2Tsv = angTsvsv(incAngSV(sv2Tsv));
    angSV2Tp  = angTsvp(incAngSV(sv2Tp));
else
    sv2sv=[]; sv2p=[]; sv2Tsv=[]; sv2Tp=[];
    angSV2SV=[]; angSV2P=[]; angSV2Tsv=[]; angSV2Tp=[];
end

% --- SH incidence: reflect-SH / transmit-SH (no mode conversion) ---
nSH = numel(incAngSH);
if nSH > 0
    Rsh = lookup(zoepAB.E_Rsh, incAngSH);
    Tsh = lookup(zoepAB.E_Tsh, incAngSH);
    Wtot = Rsh + Tsh;
    u = rand(nSH,1) .* Wtot;
    sh2sh  = u <= Rsh;
    sh2Tsh = ~sh2sh;
    angSH2SH  = angRsvsv(incAngSH(sh2sh));   % specular, same law as SV-SV
    angSH2Tsh = angTsh(incAngSH(sh2Tsh));
else
    sh2sh=[]; sh2Tsh=[];
    angSH2SH=[]; angSH2Tsh=[];
end

% --- output ---
o.pparticle  = pIdx;
o.sparticle  = sIdx;
o.svparticle = isSV;
o.shparticle = ~isSV;

o.p2p = p2p;   o.p2sv = p2sv;   o.p2Tp = p2Tp;   o.p2Tsv = p2Tsv;
o.sv2sv = sv2sv; o.sv2p = sv2p; o.sv2Tsv = sv2Tsv; o.sv2Tp = sv2Tp;
o.sh2sh = sh2sh; o.sh2Tsh = sh2Tsh;

o.angPP = angP2P;   o.angPSV = angP2SV;   o.angPTp = angP2Tp;   o.angPTsv = angP2Tsv;
o.angSVSV = angSV2SV; o.angSVP = angSV2P; o.angSVTsv = angSV2Tsv; o.angSVTp = angSV2Tp;
o.angSHSH = angSH2SH; o.angSHTsh = angSH2Tsh;

end
