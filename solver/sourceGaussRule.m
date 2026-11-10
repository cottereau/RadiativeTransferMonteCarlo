function [x,w] = sourceGaussRule(edges)
% Four-point Gauss-Legendre on each interval; positive physical weights.
z = [-0.8611363115940526;-0.3399810435848563;0.3399810435848563;0.8611363115940526];
v = [0.3478548451374539;0.6521451548625461;0.6521451548625461;0.3478548451374539];
h = diff(edges(:)')/2;
c = (edges(1:end-1)+edges(2:end))/2;
x = z*h+c(:)'; w = v*h;
x = x(:); w = w(:);
end
