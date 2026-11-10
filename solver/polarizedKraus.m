function [amp,gain] = polarizedKraus(mat,mix,a,b,n,m,e)
%POLARIZEDKRAUS Vectorized Kraus amplitudes of elasticPolarizedKernels.
% mix = polarizedKraus(mat) returns the covariance factor (lambda,mu,rho).
% a,b: 1 = P, 2 = S (incident, outgoing). n,m,e: N-by-3 incident and
% outgoing directions and incident S polarization (ignored for P).
% amp(:,:,l) is N-by-1 (P out) or N-by-3 (S out, transverse to m);
% gain = sum_l |amp_l|^2 in 1/(s sr).
if nargin==1
    cv = double(mat.coefficients_of_variation(:)); r = double(mat.correlation_coefficients(:));
    C = (cv*cv').*[1 r(1) r(2); r(1) 1 r(3); r(2) r(3) 1]; C = (C+C')/2;
    [U,D] = eig(C); amp = U*diag(sqrt(max(diag(D),0))); return
end
c = sum(n.*m,2); K = mat.vp/mat.vs;
omega = 2*pi*double(mat.Frequency);
z = omega*double(mat.CorrelationLength)./[mat.vp mat.vs];
q = sqrt(max(0,z(a)^2+z(b)^2-2*z(a)*z(b)*c));
base = (pi/2)*omega*z(1)^3;
mu = reshape(mix(2,:),1,1,3); rh = reshape(mix(3,:),1,1,3);
switch 2*(a-1)+b
    case 1 % P -> P
        A = [(1-2/K^2)*ones(size(c)), 2*c.^2/K^2, -c];
        amp = reshape(A*mix,[],1,3); pref = base;
    case 2 % P -> S
        amp = (n-c.*m).*(2*c.*mu-K*rh); pref = base*K;
    case 3 % S -> P
        amp = sum(m.*e,2).*((2/K)*c.*mu-rh); pref = base;
    case 4 % S -> S
        me = sum(m.*e,2); te = e-me.*m;
        amp = (me.*(n-c.*m)+c.*te).*mu-te.*rh; pref = (pi/2)*omega*z(2)^3;
end
amp = amp.*sqrt(pref*reshape(mat.Phi(q),[],1));
gain = reshape(sum(sum(amp.^2,2),3),[],1);
end

