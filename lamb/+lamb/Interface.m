classdef Interface
%INTERFACE  Welded plane interface z = 0 between a lower solid and an upper medium.
%
%   I = lamb.Interface(lowerSolid, upperMedium)   % 'solid-solid' or 'solid-fluid'
%   I = lamb.Interface(lowerSolid)                % 'free' (traction-free surface, Lamb)
%
%   The lower medium (z > 0) must be a solid. The boundary-condition matrix depends only on
%   the interface; the source enters through the right-hand side. Calls are dispatched to
%   lamb.SolidSolidSolver or lamb.SolidFluidSolver (fluid = [] for 'free').

    properties (SetAccess = immutable)
        lower          % lamb.Medium (solid), z > 0
        upper          % lamb.Medium or [],   z < 0
    end
    properties (Dependent)
        kind           % 'solid-solid' | 'solid-fluid' | 'free'
    end

    methods
        function obj = Interface(lower, upper)
            arguments
                lower (1,1) lamb.Medium
                upper = []
            end
            if lower.isFluid
                error('lamb:Interface:lower', 'The lower half-space must be a solid.');
            end
            if ~isempty(upper) && ~isa(upper, 'lamb.Medium')
                error('lamb:Interface:upper', 'upper must be a lamb.Medium or [].');
            end
            obj.lower = lower;
            obj.upper = upper;
        end

        function s = get.kind(obj)
            if isempty(obj.upper),     s = 'free';
            elseif obj.upper.isFluid,  s = 'solid-fluid';
            else,                      s = 'solid-solid';
            end
        end

        function name = interfaceWave(obj)
            %INTERFACEWAVE  Name of the trapped wave: 'Rayleigh' | 'Scholte' | 'Stoneley'.
            name = struct('free','Rayleigh','solid_fluid','Scholte','solid_solid','Stoneley') ...
                .(strrep(obj.kind, '-', '_'));
        end

        function kb = branchPoints(obj, omega)
            %BRANCHPOINTS  Sorted real branch points {k_a1, k_b1, k_a2, k_b2 | k_f}.
            kb = [obj.lower.wavenumber(omega,'P'), obj.lower.wavenumber(omega,'S')];
            if ~isempty(obj.upper)
                kb(end+1) = obj.upper.wavenumber(omega, 'P');
                if ~obj.upper.isFluid
                    kb(end+1) = obj.upper.wavenumber(omega, 'S');
                end
            end
            kb = sort(kb);
        end

        function M = matrix(obj, k, omega)
            M = feval([obj.solver '.matrix'], obj.lower, obj.upper, k, omega);
        end
        function S = rhs(obj, k, omega, amp)
            S = feval([obj.solver '.rhs'], obj.lower, obj.upper, k, omega, amp);
        end
        function x = solve(obj, amp, k, omega)
            %SOLVE  x = I.solve(amp, k, omega); amp from lamb.incidentAmplitudes.
            x = feval([obj.solver '.solve'], obj.lower, obj.upper, k, omega, amp);
        end
        function D = dispersion(obj, k, omega)
            D = feval([obj.solver '.dispersion'], obj.lower, obj.upper, k, omega);
        end
        function names = unknowns(obj)
            %UNKNOWNS  P-SV unknown names in matrix-column order.
            names = feval([obj.solver '.unknowns'], obj.upper);
        end
    end

    methods (Access = private)
        function s = solver(obj)
            if strcmp(obj.kind, 'solid-solid')
                s = 'lamb.SolidSolidSolver';
            else
                s = 'lamb.SolidFluidSolver';
            end
        end
    end
end
