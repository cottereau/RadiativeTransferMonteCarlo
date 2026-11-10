function [L0,L1,L2,L3] = legAssoc(n,m,u)
%% legAssoc - associated Legendre P_n^m(u) and its u-derivatives up to third order
%
% Condon-Shortley convention:  P_n^1(u) = -(1-u^2)^(1/2) * dP_n/du
%
% Implemented only for m = 0 and m = 1, which are the orders the sphere
% needs: P incidence is axisymmetric (m=0) and S incidence lives in m=1.

switch m
    case 0
        [L0,L1,L2,L3] = TMatrix.sphereFun('P',n,u);

    case 1
        if n == 0
            L0 = zeros(size(u)); L1 = L0; L2 = L0; L3 = L0; return
        end
        [~,g0,g1,g2,g3] = TMatrix.sphereFun('P',n,u);   % g = dP_n/du and its derivatives

        s   = sqrt(1-u.^2);
        s1  = -u./s;
        s2  = -1./s.^3;
        s3  = -3*u./s.^5;

        % L = -s*g, Leibniz rule
        L0 = -( s.*g0 );
        L1 = -( s1.*g0 + s.*g1 );
        L2 = -( s2.*g0 + 2*s1.*g1 + s.*g2 );
        L3 = -( s3.*g0 + 3*s2.*g1 + 3*s1.*g2 + s.*g3 );

    otherwise
        error('legAssoc is implemented only for m = 0 or 1');
end
end
