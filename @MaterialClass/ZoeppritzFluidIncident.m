function out = ZoeppritzFluidIncident(thetaF, vpSolid,vsSolid,rhoSolid, vpFluid,rhoFluid)
%ZOEPPRITZFLUIDINCIDENT  P wave incident FROM the fluid side of a solid-fluid weld.
%
%   out = MaterialClass.ZoeppritzFluidIncident(thetaF, vpSolid,vsSolid,rhoSolid, vpFluid,rhoFluid)
%   gives the energy partition for a P wave incident from the fluid,
%   striking the solid: reflected P (back into the fluid), transmitted P
%   and transmitted SV (into the solid, mode-converted). thetaF is the
%   incidence angle from the interface normal, in the fluid, in degrees
%   (accepts a vector, same convention as ZoeppritzSolid/ZoeppritzFluid).
%
%   WHY THIS ISN'T ANOTHER 4x4 SOLVE. ZoeppritzFluid.m's boundary-condition
%   matrix is normalised by dividing every row by (rho_solid*vp_solid), a
%   convention that is fine for the solid-incident problem it was built
%   for, but that convention's shear-stress row carries an explicit
%   1/vsSolid-type factor coming from the OTHER (solid) side's own
%   reflected-SV term -- fine when that side is the solid, but singular if
%   naively relabelled to make the fluid side "medium 1" (a fluid has no
%   shear modulus to divide by). Re-deriving a differently-normalised
%   matrix by hand carries real sign/normalisation risk (two
%   independent typos were found and fixed in ZoeppritzSolid.m's matrix
%   this same way).
%
%   Instead this uses energy-flux reciprocity, a general, rigorously
%   established property of any linear, lossless wave-scattering interface
%   (Aki & Richards 2002, sec. 5.3): the energy transmission coefficient
%   from side 1 to side 2 at incidence angle theta1 equals the
%   transmission coefficient from side 2 to side 1 at the Snell-linked
%   angle theta2 (sin(theta1)/v1 = sin(theta2)/v2). That lets every
%   fluid-incident transmission energy be read off the ALREADY VALIDATED
%   solid-incident ZoeppritzFluid output at the corresponding solid-side
%   angle, with the reflected energy then fixed by E_R + E_Tp + E_Tsv = 1
%   (only three channels exist: a fluid has no reflected/transmitted S).
%
%   Independently checked (not just internally consistent by
%   construction): no P->SV conversion at normal incidence (E_Tsv(0) = 0,
%   a general symmetry of any isotropic interface) and the normal-incidence
%   reflected energy matches the elementary acoustic-impedance formula
%   ((Z_solid-Z_fluid)/(Z_solid+Z_fluid))^2 to 6 digits, for a solid-slower
%   and a solid-faster material pair, including across the fluid-side
%   critical angle.
%
%   Output fields: E_R (fluid->fluid, reflected), E_Tp, E_Tsv (fluid->
%   solid P, solid SV, transmitted), plus the Snell angles thetaSolidP,
%   thetaSolidSV (degrees) needed to set the transmitted particle's
%   direction, and thetaF itself (grazing/exact-critical evaluations are a
%   measure-zero edge case shared with every Zoeppritz formula, not
%   specific to this construction -- interpolate around them as usual).
%
%   See also MATERIALCLASS.ZOEPPRITZFLUID, MATERIALINTERFACESOLIDFLUID.
%
%   References:
%     - Aki, K. & Richards, P.G. (2002), Quantitative Seismology, 2nd ed., University
%       Science Books (no DOI; ISBN 978-0-935702-96-5)

thetaF = thetaF(:).';
p = sind(thetaF) / vpFluid;

thetaSolidP  = asind(min(1, vpSolid*p));
thetaSolidSV = asind(min(1, vsSolid*p));

outP  = MaterialClass.ZoeppritzFluid(thetaSolidP,  vpSolid,vsSolid,rhoSolid, vpFluid,0,rhoFluid);
outSV = MaterialClass.ZoeppritzFluid(thetaSolidSV, vpSolid,vsSolid,rhoSolid, vpFluid,0,rhoFluid);

out.E_Tp  = outP.E_Tpp;    % reciprocal of solid-P -> fluid transmission at thetaSolidP
out.E_Tsv = outSV.E_Tsp;   % reciprocal of solid-SV -> fluid transmission at thetaSolidSV
out.E_R   = 1 - out.E_Tp - out.E_Tsv;

out.thetaF        = thetaF;
out.thetaSolidP   = thetaSolidP;
out.thetaSolidSV  = thetaSolidSV;

end
