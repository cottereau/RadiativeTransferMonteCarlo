% Runs every reflection/transmission-coefficient example in this folder
% and prints a one-line PASS/FAIL summary for each.
%
% Every example now runs TWO phases at the same material pair and the same
% incidence angles:
%
%   [canal]  the interface operator alone. MaterialInterfaceTwoSolid /
%            MaterialInterfaceSolidFluid is called once on a synthetic
%            particle set sitting exactly on the interface, and the fraction
%            sent into each channel is compared with the Zoeppritz energy
%            coefficients.
%
%   [RTE ]   the transport end to end. A finite plane-wave packet is launched
%            away from the interface inside medium 1, propagateParticleMultiMaterial
%            carries it across, and the energy is measured afterwards from
%            where the particles actually ended up -- which side they are on,
%            which mode they carry, and how their polarization projects onto
%            the SV/SH basis. Each channel is reported as the level it SHOULD
%            settle at against the level it ACTUALLY reaches, and the energy
%            history is plotted with those levels drawn over it.
%
% The second phase is what catches anything only the transport can break: a
% lost or duplicated particle, a matIdx left stale after a crossing, a
% polarization that does not survive the interface. The first phase cannot
% see any of those, because it never propagates anything.
%
%   solid-solid:  TwoSolid_Pincidence, TwoSolid_SVincidence, TwoSolid_SHincidence
%   solid-fluid:  SolidFluid_Pincidence, SolidFluid_SVincidence (solid incident)
%                 FluidSolid_Pincidence (fluid incident, via reciprocity)
%
% See also MATERIALINTERFACETWOSOLID, MATERIALINTERFACESOLIDFLUID.

clc
tests = {
    @TwoSolid_Pincidence
    @TwoSolid_SVincidence
    @TwoSolid_SHincidence
    @SolidFluid_Pincidence
    @SolidFluid_SVincidence
    @FluidSolid_Pincidence
};

results = false(size(tests));
for k = 1:numel(tests)
    results(k) = tests{k}();
end

fprintf('\n========================================\n');
fprintf('InterfaceCoefficients summary:\n');
for k = 1:numel(tests)
    name = func2str(tests{k});
    if results(k), tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  %-24s [%s]\n', name, tag);
end
if all(results)
    fprintf('ALL EXAMPLES PASSED\n');
else
    fprintf('SOME EXAMPLES FAILED -- see the per-channel output above\n');
end
fprintf('========================================\n');
