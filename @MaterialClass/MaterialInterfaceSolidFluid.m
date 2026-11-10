function o = MaterialInterfaceSolidFluid(matA, fromFluid, zoepAB, P, n, ind)
%MATERIALINTERFACESOLIDFLUID  Reflect/transmit particles at a solid-fluid weld.
%
%   o = MaterialInterfaceSolidFluid(matA, fromFluid, zoepAB, P, n, ind)
%   mirrors MaterialInterfaceTwoSolid but for an interface where one side
%   is a fluid (no shear modulus, no S waves at all): matA is the medium
%   the incident particles ind are currently in, fromFluid says whether
%   that medium is the fluid or the solid side, and zoepAB is the
%   precomputed table for THIS side -- from MaterialClass.ZoeppritzFluid
%   (solid-incident) or MaterialClass.ZoeppritzFluidIncident (fluid-
%   incident), built once per material pair by the caller.
%
%   Every particle currently in the fluid is a P wave (P.p is true): a
%   fluid supports no S wave to be, so this is a domain invariant rather
%   than something checked here. Outcomes: from the fluid, reflect-P or
%   transmit into the solid as P or (mode-converted) SV. From the solid,
%   P and SV incidence each give reflect-same / reflect-converted /
%   transmit-P-into-fluid (a fluid has no transmitted-S channel to speak
%   of); SH incidence always reflects (a fluid cannot support any shear
%   traction, so SH reflects with R = 1, same as at a free surface).
%
%   See also MATERIALINTERFACETWOSOLID, ZOEPPRITZFLUID, ZOEPPRITZFLUIDINCIDENT.

if fromFluid
    o = fluidSideDecision(zoepAB, P, n, ind);
else
    o = solidSideDecision(matA, zoepAB, P, n, ind);
end

end

% -------------------------------------------------------------------------
function o = fluidSideDecision(zoepAB, P, n, ind)
% every particle in ind is a P wave (fluid has no S)
partP_dir = P.dir(ind,:);
nP = size(partP_dir,1);
o.pparticle  = ind;
o.sparticle  = false(P.N,1);
o.svparticle = false(0,1);
o.shparticle = false(0,1);

if nP == 0
    o.p2p = []; o.p2Tp = []; o.p2Tsv = [];
    o.angPP = []; o.angPTp = []; o.angPTsv = [];
    return
end

nRep = n(ind,:);
incAng = acosd(min(1,abs(dot(partP_dir, nRep, 2))));

lookup = @(data, ang) interp1(zoepAB.thetaF, data, ang, 'linear', 'extrap');
lookupAng = @(data, ang) interp1(zoepAB.thetaF, data, ang, 'linear', 'extrap');

Rf  = lookup(zoepAB.E_R,   incAng);
Tp  = lookup(zoepAB.E_Tp,  incAng);
Tsv = lookup(zoepAB.E_Tsv, incAng);
Wtot = Rf + Tp + Tsv;
u = rand(nP,1) .* Wtot;
c1 = Rf; c2 = c1+Tp;
p2p   = u <= c1;
p2Tp  = u > c1 & u <= c2;
p2Tsv = u > c2;

o.p2p   = p2p;
o.p2Tp  = p2Tp;
o.p2Tsv = p2Tsv;

o.angPP   = incAng(p2p);                                  % specular
o.angPTp  = lookupAng(zoepAB.thetaSolidP,  incAng(p2Tp));
o.angPTsv = lookupAng(zoepAB.thetaSolidSV, incAng(p2Tsv));
end

% -------------------------------------------------------------------------
function o = solidSideDecision(matA, zoepAB, P, n, ind)
% zoepAB has fields from both incidence tables, built by the caller as
% zoepAB.P (MaterialClass.ZoeppritzFluid, P incidence) and zoepAB.SV
% (same, SV incidence); both share the same theta grid, zoepAB.j1_deg.
j1_deg = zoepAB.j1_deg;
lookup = @(data, ang) interp1(j1_deg, data, ang, 'linear', 'extrap');

vp = matA.vp; vs = matA.vs;
angRpp  = @(th) th;
angRpsv = @(th) asind(min(1,(sind(th)/vp)*vs));
angTp   = @(th, mode) th; % filled per-branch below (P or SV incidence has its own transmitted-P Snell map)

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

% --- P incidence: reflect-P / reflect-SV / transmit-P(into fluid) ---
if nP > 0
    Rpp  = lookup(zoepAB.P.E_Rpp,  incAngP);
    Rpsv = lookup(zoepAB.P.E_Rpsv, incAngP);
    Tpp  = lookup(zoepAB.P.E_Tpp,  incAngP);
    Wtot = Rpp + Rpsv + Tpp;
    u = rand(nP,1) .* Wtot;
    c1 = Rpp; c2 = c1+Rpsv;
    p2p  = u <= c1;
    p2sv = u > c1 & u <= c2;
    p2Tp = u > c2;

    angP2P  = angRpp(incAngP(p2p));
    angP2SV = angRpsv(incAngP(p2sv));
    angP2Tp = asind(min(1,(sind(incAngP(p2Tp))/vp)*zoepAB.vpFluid));
else
    p2p=[]; p2sv=[]; p2Tp=[];
    angP2P=[]; angP2SV=[]; angP2Tp=[];
end

% --- SV incidence: reflect-SV / reflect-P / transmit-P(into fluid) ---
nSV = numel(incAngSV);
if nSV > 0
    Rsvsv = lookup(zoepAB.SV.E_Rsvsv, incAngSV);
    Rsp   = lookup(zoepAB.SV.E_Rsp,   incAngSV);
    Tsp   = lookup(zoepAB.SV.E_Tsp,   incAngSV);
    Wtot = Rsvsv + Rsp + Tsp;
    u = rand(nSV,1) .* Wtot;
    c1 = Rsvsv; c2 = c1+Rsp;
    sv2sv = u <= c1;
    sv2p  = u > c1 & u <= c2;
    sv2Tp = u > c2;

    angSV2SV = angRpp(incAngSV(sv2sv));                               % specular
    angSV2P  = asind(min(1,(sind(incAngSV(sv2p))/vs)*vp));
    angSV2Tp = asind(min(1,(sind(incAngSV(sv2Tp))/vs)*zoepAB.vpFluid));
else
    sv2sv=[]; sv2p=[]; sv2Tp=[];
    angSV2SV=[]; angSV2P=[]; angSV2Tp=[];
end

% --- SH incidence: a fluid cannot carry any shear traction, so SH always
% reflects (R = 1), exactly as at a free surface.
nSH = numel(incAngSH);
if nSH > 0
    sh2sh = true(nSH,1);
    angSH2SH = angRpp(incAngSH);
else
    sh2sh = []; angSH2SH = [];
end

o.pparticle  = pIdx;
o.sparticle  = sIdx;
o.svparticle = isSV;
o.shparticle = ~isSV;

o.p2p = p2p;   o.p2sv = p2sv;   o.p2Tp = p2Tp;
o.sv2sv = sv2sv; o.sv2p = sv2p; o.sv2Tp = sv2Tp;
o.sh2sh = sh2sh;

o.angPP = angP2P;   o.angPSV = angP2SV;   o.angPTp = angP2Tp;
o.angSVSV = angSV2SV; o.angSVP = angSV2P; o.angSVTp = angSV2Tp;
o.angSHSH = angSH2SH;
end
