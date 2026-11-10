function D = sourceMomentDeviator(M,scale)
% Scaled deviator retaining small differences between near-equal diagonals.
D = M/scale;
d12 = difference(M(1,1),M(2,2),scale);
d13 = difference(M(1,1),M(3,3),scale);
d23 = difference(M(2,2),M(3,3),scale);
D(1,1) = (d12+d13)/3;
D(2,2) = (-d12+d23)/3;
D(3,3) = (-d13-d23)/3;
end

function d = difference(a,b,scale)
if sign(a)==sign(b), d = (a-b)/scale;
else, d = a/scale-b/scale; end
end
