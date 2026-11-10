% References:
%   - Aki, K. & Richards, P.G. (2002), Quantitative Seismology, 2nd ed., University Science
%     Books (no DOI; ISBN 978-0-935702-96-5)

classdef PointSource
%POINTSOURCE  Point force of arbitrary direction or arbitrary moment tensor at (0,0,zs).
%
%   s = lamb.PointSource('force',  F, zs)     % F  3-by-1 [N], axes (x, y, z-down)
%   s = lamb.PointSource('moment', M, zs)     % M  3-by-3 [N m] (symmetric = seismic)
%   s = lamb.PointSource('general', s9, zs)   % force AND moment together (e.g. the Born
%                                             % equivalent source of a scatterer), s9 the
%                                             % 9-vector of lamb.rte.sourceVector; may be complex
%   s = lamb.PointSource(..., 'Wavelet', @(omega) ...)   % spectrum s(omega); must accept
%                                                        % complex omega (default 1)
%   Helpers: lamb.PointSource.explosion(M0, zs), .doubleCouple(M0, strike, dip, rake, zs),
%            lamb.PointSource.ricker(fp, t0)  (wavelet handle)
%
%   zs >= 0 : source in the lower solid (zs = 0 means just below the interface).
%   zs <  0 : source in the upper medium (fluid only, Case 3).

    properties (SetAccess = immutable)
        kind    char = 'force'
        F       (3,1) double = [0; 0; 1]
        M       (3,3) double = zeros(3)
        zs      (1,1) double = 0
        wavelet function_handle = @(w) ones(size(w))
    end

    methods
        function obj = PointSource(kind, FM, zs, opts)
            arguments
                kind {mustBeMember(kind, {'force','moment','general'})}
                FM double
                zs (1,1) double
                opts.Wavelet function_handle = @(w) ones(size(w))
            end
            obj.kind = char(kind);
            switch obj.kind
                case 'force',  obj.F = FM(:);
                case 'moment', obj.M = FM;
                otherwise,     B = lamb.rte.sourceVector(FM(:));  obj.F = B.F;  obj.M = B.M;
            end
            obj.zs = zs;
            obj.wavelet = opts.Wavelet;
        end

        function [Pi, Sig, H] = directAmplitudes(obj, med, k, phi, omega, epsSign)
            %DIRECTAMPLITUDES  Whole-space amplitudes of e^{-nu|z-zs|} (eq:srcF, eq:srcM, eq:PSH).
            %   [Pi, Sig, H] = src.directAmplitudes(med, k, phi, omega, epsSign)
            %   k N-by-1, phi 1-by-Nphi, epsSign = sgn(z - zs) = +1 (below) or -1 (above).
            %   Pi : P potential, Sig : SV potential, H : SH displacement (N-by-Nphi).
            %   In a fluid only Pi is non-zero.
            k = k(:);  phi = phi(:).';
            ikx = 1i*k.*cos(phi);  iky = 1i*k.*sin(phi);
            rw  = med.rho*omega.^2;
            na  = med.verticalWavenumber(k, omega, 'P');
            kap = {ikx, iky, -epsSign*na.*ones(size(phi))};           % kappa_alpha
            switch obj.kind                                            % linear: general = sum
                case 'force',  cP = -obj.dotF(kap)./rw;
                case 'moment', cP = obj.quadM(kap)./rw;
                otherwise,     cP = (obj.quadM(kap) - obj.dotF(kap))./rw;
            end
            Pi = cP./(2*na);
            if med.isFluid
                Sig = zeros(size(Pi));  H = Sig;
                return
            end
            nb  = med.verticalWavenumber(k, omega, 'S');
            kb2 = med.wavenumber(omega, 'S').^2;
            kbv = {ikx, iky, -epsSign*nb.*ones(size(phi))};           % kappa_beta
            tx = -sin(phi);  ty = cos(phi);                            % t-hat (kappa.t = 0)
            wz = 0;  wt = 0;
            if ~strcmp(obj.kind, 'moment')
                wz = (kb2.*obj.F(3) + kbv{3}.*obj.dotF(kbv))./rw;
                wt = kb2.*(obj.F(1)*tx + obj.F(2)*ty)./rw;
            end
            if ~strcmp(obj.kind, 'force')
                Mk = obj.Mkappa(kbv);
                wz = wz - (kb2.*Mk{3} + kbv{3}.*obj.quadM(kbv))./rw;
                wt = wt - kb2.*(tx.*Mk{1} + ty.*Mk{2})./rw;
            end
            Sig = wz./(2i*k.*nb);
            H   = wt./(2*nb);
        end
    end

    methods (Access = private)
        function v = dotF(obj, kap)          % kappa . F
            v = kap{1}*obj.F(1) + kap{2}*obj.F(2) + kap{3}*obj.F(3);
        end
        function Mk = Mkappa(obj, kap)       % (M kappa)_p
            Mk = cell(1,3);
            for p = 1:3
                Mk{p} = obj.M(p,1)*kap{1} + obj.M(p,2)*kap{2} + obj.M(p,3)*kap{3};
            end
        end
        function v = quadM(obj, kap)         % kappa^T M kappa (bilinear)
            Mk = obj.Mkappa(kap);
            v = kap{1}.*Mk{1} + kap{2}.*Mk{2} + kap{3}.*Mk{3};
        end
    end

    methods (Static)
        function s = explosion(M0, zs, varargin)
            s = lamb.PointSource('moment', M0*eye(3), zs, varargin{:});
        end

        function s = doubleCouple(M0, strike, dip, rake, zs, varargin)
            %DOUBLECOUPLE  Aki & Richards (Box 4.4), x = north, y = east, z = down; angles in deg.
            f = deg2rad(strike); d = deg2rad(dip); l = deg2rad(rake);
            Mxx = -(sin(d)*cos(l)*sin(2*f) + sin(2*d)*sin(l)*sin(f)^2);
            Mxy =   sin(d)*cos(l)*cos(2*f) + 0.5*sin(2*d)*sin(l)*sin(2*f);
            Mxz = -(cos(d)*cos(l)*cos(f)   + cos(2*d)*sin(l)*sin(f));
            Myy =   sin(d)*cos(l)*sin(2*f) - sin(2*d)*sin(l)*cos(f)^2;
            Myz = -(cos(d)*cos(l)*sin(f)   - cos(2*d)*sin(l)*cos(f));
            Mzz =   sin(2*d)*sin(l);
            M = M0*[Mxx Mxy Mxz; Mxy Myy Myz; Mxz Myz Mzz];
            s = lamb.PointSource('moment', M, zs, varargin{:});
        end

        function W = ricker(fp, t0)
            %RICKER  Spectrum of r(t) = (1 - 2 pi^2 fp^2 tau^2) exp(-pi^2 fp^2 tau^2), tau = t - t0,
            %   with the forward transform int r(t) e^{+i omega t} dt (e^{-i omega t} convention).
            wp = 2*pi*fp;
            W = @(w) (4*sqrt(pi)/wp^3)*w.^2.*exp(-(w/wp).^2).*exp(1i*w*t0);
        end
    end
end
