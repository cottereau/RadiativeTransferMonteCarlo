function R = rotationMatrix(phi)
%ROTATIONMATRIX  9-by-9 action on s of a rotation by phi about the z axis.
%   R = lamb.rte.rotationMatrix(phi):  s' = R*s  is the source rotated by +phi
%   (F' = Rz F, M' = Rz M Rz.').
%
%   The medium is invariant about z, so the amplitude radiated at azimuth phi by s equals the
%   amplitude radiated at azimuth 0 by the source rotated by -phi:
%       D(phi; s) = g(0).' * rotationMatrix(-phi) * s.
%   Tables are therefore stored at phi = 0 only.
c = cos(phi);  s = sin(phi);
Rz = [c -s 0; s c 0; 0 0 1];
R = zeros(9);
for i = 1:9
    e = zeros(9, 1);  e(i) = 1;
    u = lamb.rte.sourceVector(e);
    R(:, i) = lamb.rte.sourceVector(Rz*u.F, Rz*u.M*Rz.');
end
end
