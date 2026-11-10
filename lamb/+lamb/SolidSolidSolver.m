classdef SolidSolidSolver
%SOLIDSOLIDSOLVER  f-k boundary-condition system for a welded solid-solid interface
%   (Case 1, eq:case1). Solid 1 (z > 0) contains the source; solid 2 occupies z < 0.
%
%   P-SV : M1 [A1; B1; A2; B2] = S1
%          rows    u_x, u_z, sigma_zz/mu1, sigma_xz/mu1
%          columns A1 (P1 refl), B1 (SV1 refl), A2 (P2 trans), B2 (SV2 trans)
%   SH   : closed form of the 2x2 system,
%          C1 = (mu1 nu_b1 - mu2 nu_b2)/(mu1 nu_b1 + mu2 nu_b2) C0,
%          C2 = 2 mu1 nu_b1 /(mu1 nu_b1 + mu2 nu_b2) C0.
%          The SH determinant never vanishes, so there is no SH interface wave.
%
%   Vectorisation: as lamb.SolidFluidSolver (k N-by-1, amplitudes N-by-Nphi, pagemldivide).

    methods (Static)
        function M = matrix(s1, s2, k, omega)
            %MATRIX  4-by-4-by-N pages.
            q  = lamb.SolidSolidSolver.coeffs(s1, s2, k, omega);
            N  = numel(q.k);
            pg = @(v) reshape(v.*ones(N,1), 1, 1, N);
            ik = 1i*q.k;  m = q.m;
            M = [pg(ik),            pg(q.nb1),         pg(-ik),               pg(q.nb2)
                 pg(-q.na1),        pg(ik),            pg(-q.na2),            pg(-ik)
                 pg(q.G1),          pg(-2*ik.*q.nb1),  pg(-m*q.G2),           pg(-2*m*ik.*q.nb2)
                 pg(-2*ik.*q.na1),  pg(-q.G1),         pg(-2*m*ik.*q.na2),    pg(m*q.G2)];
        end

        function S = rhs(s1, s2, k, omega, amp)
            %RHS  4-by-Nphi-by-N right-hand sides from the up-going source field (P0, S0).
            if isfield(amp, 'Pf0')
                error('lamb:SolidSolidSolver:source', 'The source must be in solid 1 (zs >= 0).');
            end
            q  = lamb.SolidSolidSolver.coeffs(s1, s2, k, omega);
            pg = @(v) permute(v, [3 2 1]);
            ik = 1i*q.k;
            P0 = amp.P0;  S0 = amp.S0;
            S = -[pg(ik.*P0 - q.nb1.*S0)
                  pg(q.na1.*P0 + ik.*S0)
                  pg(q.G1.*P0 + 2*ik.*q.nb1.*S0)
                  pg(2*ik.*q.na1.*P0 - q.G1.*S0)];
        end

        function x = solve(s1, s2, k, omega, amp)
            %SOLVE  Scattered amplitudes, each N-by-Nphi: .A1 .B1 .A2 .B2 .C1 .C2
            M = lamb.SolidSolidSolver.matrix(s1, s2, k, omega);
            X = pagemldivide(M, lamb.SolidSolidSolver.rhs(s1, s2, k, omega, amp));
            names = lamb.SolidSolidSolver.unknowns();
            for i = 1:4
                x.(names{i}) = permute(X(i,:,:), [3 2 1]);
            end
            q = lamb.SolidSolidSolver.coeffs(s1, s2, k, omega);
            den  = q.nb1 + q.m*q.nb2;
            x.C1 = (q.nb1 - q.m*q.nb2)./den .* amp.C0;
            x.C2 = 2*q.nb1./den .* amp.C0;
        end

        function D = dispersion(s1, s2, k, omega)
            %DISPERSION  det(M1) (Stoneley equation); real for k beyond all branch points.
            M = lamb.SolidSolidSolver.matrix(s1, s2, k, omega);
            D = zeros(size(M,3), 1);
            for i = 1:numel(D)
                D(i) = det(M(:,:,i));
            end
        end

        function names = unknowns(~)
            names = {'A1','B1','A2','B2'};
        end
    end

    methods (Static, Access = private)
        function q = coeffs(s1, s2, k, omega)
            q.k   = k(:);
            q.na1 = s1.verticalWavenumber(q.k, omega, 'P');
            q.nb1 = s1.verticalWavenumber(q.k, omega, 'S');
            q.na2 = s2.verticalWavenumber(q.k, omega, 'P');
            q.nb2 = s2.verticalWavenumber(q.k, omega, 'S');
            q.G1  = 2*q.k.^2 - s1.wavenumber(omega, 'S').^2;
            q.G2  = 2*q.k.^2 - s2.wavenumber(omega, 'S').^2;
            q.m   = s2.mu/s1.mu;
        end
    end
end
