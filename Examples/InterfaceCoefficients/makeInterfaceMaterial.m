function mat = makeInterfaceMaterial(vp, vs, rho)
%MAKEINTERFACEMATERIAL  Homogeneous 3D material with no volume scattering.
%
%   mat = makeInterfaceMaterial(vp, vs, rho) builds the minimal MaterialClass
%   propagateParticleMultiMaterial needs, with meanFreeTime = Inf so that the
%   only thing that can change a particle's mode or direction is the
%   interface. Pass vs = 0 for a fluid.
%
%   This is the two-material counterpart of the local helper in
%   Examples/BoundaryReflectionEnergyTests.m, kept shared because all six
%   InterfaceCoefficients examples need exactly the same thing.

mat = MaterialClass();
mat.d = 3;
mat.acoustics = false;
mat.vp = vp;
mat.vs = vs;
mat.rho = rho;
mat.timeSteps = 0;

% No volume scattering: the interface is the only physics under test.
mat.meanFreeTime = [Inf; Inf];
mat.P2P = 1;
mat.S2S = 1;

zeroInv = @(u) zeros(size(u));
mat.invcdf = cell(2,2);
mat.invcdf{1,1} = zeroInv;
mat.invcdf{1,2} = zeroInv;
mat.invcdf{2,1} = zeroInv;
mat.invcdf{2,2} = zeroInv;
end
