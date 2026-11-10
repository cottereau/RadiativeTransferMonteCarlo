classdef SolidFluidSolver
%SOLIDFLUIDSOLVER  f-k boundary-condition systems: solid (z > 0) below a fluid (z < 0),
%   or below vacuum (traction-free surface, fluid = []).
%
%   Cases (LambFK_Formulation.tex):
%     Case 2  source in the solid   M2 [A1; B1; Af] = S2          (eq:case2)
%     Case 3  source in the fluid   M2 [A1; B1; Af] = S3          (eq:case3)  same M2
%     Free    fluid = []            M_FS [A1; B1]   = S2(2:3)     (Lamb / Rayleigh)
%     SH      sigma_zT = 0  ->  C1 = C0 (total reflection); a fluid source radiates no SH.
%
%   Rows    : u_z, sigma_zz/mu, sigma_xz/mu  (stress rows scaled by 1/mu of the solid)
%   Columns : A1 (P, e^{-nu_a z}), B1 (SV, e^{-nu_b z}), Af (fluid P, e^{+nu_f z})
%
%   Vectorisation: k is N-by-1, omega is a scalar or N-by-1 (real or complex), the source
%   amplitudes are N-by-Nphi. M depends on (k, omega) only, so every page is factorised once
%   for all Nphi right-hand sides (pagemldivide).

    methods (Static)
        function M = matrix(solid, fluid, k, omega)
            %MATRIX  n-by-n-by-N pages, n = 3 (fluid) or 2 (fluid = []).
            q  = lamb.SolidFluidSolver.coeffs(solid, fluid, k, omega);
            N  = numel(q.k);
            pg = @(v) reshape(v.*ones(N,1), 1, 1, N);
            ik = 1i*q.k;
            if isempty(fluid)
                M = [pg(q.G),            pg(-2*ik.*q.nb)
                     pg(-2*ik.*q.na),    pg(-q.G)];
            else
                M = [pg(-q.na),          pg(ik),            pg(-q.nf)
                     pg(q.G),            pg(-2*ik.*q.nb),   pg(q.rw)
                     pg(-2*ik.*q.na),    pg(-q.G),          pg(0)];
            end
        end

        function S = rhs(solid, fluid, k, omega, amp)
            %RHS  n-by-Nphi-by-N right-hand sides.
            %   amp.P0, amp.S0 (source in solid, up-going at the interface)  -> Case 2 / free
            %   amp.Pf0        (source in fluid, down-going at the interface) -> Case 3
            q  = lamb.SolidFluidSolver.coeffs(solid, fluid, k, omega);
            pg = @(v) permute(v, [3 2 1]);          % N-by-Nphi -> 1-by-Nphi-by-N
            ik = 1i*q.k;
            if isfield(amp, 'Pf0')
                if isempty(fluid)
                    error('lamb:SolidFluidSolver:noFluid', 'A source above a free surface is not defined.');
                end
                S = [pg(-q.nf.*amp.Pf0); pg(-q.rw.*amp.Pf0); pg(zeros(size(amp.Pf0)))];
            else
                s2 = -(q.G.*amp.P0 + 2*ik.*q.nb.*amp.S0);
                s3 = -(2*ik.*q.na.*amp.P0 - q.G.*amp.S0);
                if isempty(fluid)
                    S = [pg(s2); pg(s3)];
                else
                    s1 = -(q.na.*amp.P0 + ik.*amp.S0);
                    S = [pg(s1); pg(s2); pg(s3)];
                end
            end
        end

        function x = solve(solid, fluid, k, omega, amp)
            %SOLVE  Scattered amplitudes, each N-by-Nphi: .A1 .B1 [.Af] .C1
            M = lamb.SolidFluidSolver.matrix(solid, fluid, k, omega);
            X = pagemldivide(M, lamb.SolidFluidSolver.rhs(solid, fluid, k, omega, amp));
            names = lamb.SolidFluidSolver.unknowns(fluid);
            for i = 1:numel(names)
                x.(names{i}) = permute(X(i,:,:), [3 2 1]);
            end
            if isfield(amp, 'C0')
                x.C1 = amp.C0;                      % SH: total reflection at sigma_zT = 0
            else
                x.C1 = zeros(size(x.A1));           % fluid source: no SH
            end
        end

        function D = dispersion(solid, fluid, k, omega)
            %DISPERSION  Dispersion function (complex on the real axis inside the windows; real
            %   for k beyond all branch points).
            %   free : R(k) = Gamma^2 - 4 k^2 nu_a nu_b                       (Rayleigh)
            %   fluid: nu_f R(k) + (rho_f/rho) k_b^4 nu_a  = det(M2)/mu^2     (Scholte, eq:scholte)
            q = lamb.SolidFluidSolver.coeffs(solid, fluid, k, omega);
            R = q.G.^2 - 4*q.k.^2.*q.na.*q.nb;
            if isempty(fluid)
                D = R;
            else
                D = q.nf.*R + (fluid.rho/solid.rho)*q.kb2.^2.*q.na;
            end
        end

        function names = unknowns(fluid)
            if isempty(fluid), names = {'A1','B1'}; else, names = {'A1','B1','Af'}; end
        end
    end

    methods (Static, Access = private)
        function q = coeffs(solid, fluid, k, omega)
            q.k   = k(:);
            q.na  = solid.verticalWavenumber(q.k, omega, 'P');
            q.nb  = solid.verticalWavenumber(q.k, omega, 'S');
            q.kb2 = solid.wavenumber(omega, 'S').^2;
            q.G   = 2*q.k.^2 - q.kb2;                          % Gamma = 2k^2 - k_beta^2
            if ~isempty(fluid)
                q.nf = fluid.verticalWavenumber(q.k, omega, 'P');
                q.rw = fluid.rho*omega.^2/solid.mu;            % rho_f omega^2 / mu
            end
        end
    end
end
