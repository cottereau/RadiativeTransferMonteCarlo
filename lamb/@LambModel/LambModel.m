classdef LambModel
%LAMBMODEL  Lamb (half-space Green's function) Monte Carlo model of one planar interface.
%   M = LambModel(iface, omega)
%   M = LambModel(iface, omega, 'dsrc', d, 'zetaMin', 1e-3, 'NPhiBins', 360)
%   M = LambModel(S)        S: model struct saved by older code (load(file).M)
%
%   iface  lamb.Interface (lower solid at z > 0; upper solid, fluid or vacuum at z < 0)
%   omega  angular frequency of the MC run
%   dsrc   distance below which a source uses the Lamb emission K1 (default 3 lambda_S of
%          the lower solid). Beyond it: whole-space emission + Zoeppritz crossing.
%   zetaMin minimum k_beta1*|z| of MOMENT sources
%
%   Value class: assigning M to another variable copies it, as the struct it replaces.
%   An EMPTY property means "not set" (e.g. het empty = no scattering model yet).
%
%   Methods (same names and signatures as the former functions):
%     lambScatteringModel   add heterogeneity and scattering tables
%     initializeParticlesLamb, propagateParticleLamb, propagateParticleCylinder, observeLamb
%     lambExportCpp         HDF5 export for the C++ engine
%     lambCancel, lambEmitLocal, lambLocalForms, lambLocalTables, lambPowerForm,
%     lambRandomGram, lambRandomForms, lambStationarity   (local slab / random-medium kernels)
%
%   See also LAMBRTEMODEL, LAMB.INTERFACE.

    properties
        % --- built by the constructor ---
        iface                 % lamb.Interface
        omega                 % angular frequency of the MC run
        kb1                   % k_beta1 = omega/beta of the lower solid
        dsrc                  % K1 emission distance
        zetaMin               % minimum k_beta1*|z| of moment sources
        T                     % Lamb radiation tables (lamb.rte.radiationTables)
        Tws                   % whole-space tables of the lower solid
        ifaceFlip             % solid-solid: the interface seen from the upper side
        Tflip                 % ... and its tables
        dslab                 % slab half-thickness of the local transport
        phiBins               % azimuth bins
        Rneg                  % 9-by-9-by-numel(phiBins), rotationMatrix(-phi)
        chan                  % channels: .name .side (1: z>0, 2: z<0) .mode (1 P, 2 S) .c
        iw                    % interface waves: .c speed
        leaky                 % leaky waves: .c speed, .rate per unit length, .sinf fluid angle
        % --- added by lambScatteringModel ---
        het                   % {hetLower, hetUpper}; empty = no scattering
        bulk                  % lambBulkScatter totals per side
        iwT                   % interface-wave scattering tables
        loc                   % local slab tables (lambLocalTables)
        cancel                % lambCancel options (signed Monte Carlo)
        zMirror               % optional specular walls [zTop zBottom]
        % --- cylinder engine (propagateParticleCylinder) ---
        geom                  % cylinder geometry (.R, .h)
        edge                  % edge coefficient tables
        noIW                  % true: P/S reflections only, no edge capture
    end

    methods
        function M = LambModel(varargin)
            if nargin == 0
                return
            end
            if nargin == 1 && (isstruct(varargin{1}) || isa(varargin{1}, 'LambModel'))
                S = varargin{1};
                if isa(S, 'LambModel'), M = S; return, end
                f = fieldnames(S);
                bad = setdiff(f, properties(M));
                if ~isempty(bad)
                    error('LambModel:unknownField', 'Unknown model field(s): %s', strjoin(bad, ', '));
                end
                for i = 1:numel(f)
                    M.(f{i}) = S.(f{i});
                end
                return
            end
            M = build(M, varargin{:});
        end
    end

    methods (Access = private)
        function M = build(M, iface, omega, opts)
            arguments
                M
                iface (1,1) lamb.Interface
                omega (1,1) double {mustBePositive}
                opts.dsrc     (1,1) double = NaN
                opts.zetaMin  (1,1) double = 1e-3
                opts.NPhiBins (1,1) double = 360
            end
            M.iface = iface;  M.omega = omega;
            M.kb1 = omega/iface.lower.beta;
            M.dsrc = opts.dsrc;
            if isnan(M.dsrc), M.dsrc = 3*2*pi/M.kb1; end
            M.zetaMin = opts.zetaMin;
            M.T   = lamb.rte.radiationTables(iface);
            M.Tws = lamb.rte.radiationTables(lamb.Interface(iface.lower, iface.lower));
            M.ifaceFlip = [];  M.Tflip = [];
            if strcmp(iface.kind, 'solid-solid')
                M.ifaceFlip = lamb.Interface(iface.upper, iface.lower);
                M.Tflip = lamb.rte.radiationTables(M.ifaceFlip);  % sources/scatterers in the upper solid
            end
            M.dslab = M.dsrc;                                    % slab half-thickness (lambScatteringModel)

            M.phiBins = 2*pi*((1:opts.NPhiBins) - 0.5)/opts.NPhiBins;
            M.Rneg = zeros(9, 9, opts.NPhiBins);
            for m = 1:opts.NPhiBins
                M.Rneg(:, :, m) = lamb.rte.rotationMatrix(-M.phiBins(m));
            end

            L = iface.lower;  U = iface.upper;
            ch = {'P1', 1, 1, L.alpha;  'SV1', 1, 2, L.beta;  'SH1', 1, 2, L.beta};
            if ~isempty(U)
                if U.isFluid
                    ch(end+1, :) = {'Pf', 2, 1, U.alpha};
                else
                    ch(end+1:end+3, :) = {'P2', 2, 1, U.alpha;  'SV2', 2, 2, U.beta;  'SH2', 2, 2, U.beta};
                end
            end
            M.chan = cell2struct(ch, {'name', 'side', 'mode', 'c'}, 2).';

            w0 = M.T.omega0;                                     % tables are at omega0 (k_beta1 = 1)
            M.iw = struct('c', num2cell(w0./[M.T.poles.k]));      % non-dispersive: c = omega0/k_p
            M.leaky = struct('c', {}, 'rate', {}, 'sinf', {});
            for j = 1:numel(M.T.leaky)
                Lk = M.T.leaky(j);
                M.leaky(j).c    = w0/Lk.a;
                M.leaky(j).rate = 2*Lk.b*omega/w0;                % energy decay per unit length, 2 Im k_L
                M.leaky(j).sinf = Lk.a*U.alpha/w0;                % sin(theta_f), theta_f from the normal
            end
        end
    end
end
