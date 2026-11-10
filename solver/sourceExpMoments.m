function J = sourceExpMoments(z,degree)
% J(:,k+1) = integral_0^1 u^k exp(i*z*u) du. Stable at z=0.
z = z(:); J = complex(zeros(numel(z),degree+1));
small = abs(z)<=2;
for k = 0:degree
    term = ones(nnz(small),1); value = term/(k+1);
    for m = 1:32
        term = term.*(1i*z(small))/m;
        value = value+term/(m+k+1);
    end
    J(small,k+1) = value;
end
a = 1i*z(~small);
J(~small,1) = expm1(a)./a;
for k = 1:degree
    J(~small,k+1) = exp(a)./a-k*J(~small,k)./a;
end
end
