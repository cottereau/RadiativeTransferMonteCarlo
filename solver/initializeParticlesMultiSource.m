function P = initializeParticlesMultiSource(sources, Ntotal, d, acoustics)
%INITIALIZEPARTICLESMULTISOURCE  Initialize particles from multiple weighted sources.
%
%   P = initializeParticlesMultiSource(sources, Ntotal, d, acoustics)
%   generates a particle population distributed across a heterogeneous collection
%   of sources, each carrying an energy fraction specified in its .fraction field.
%
%   INPUTS:
%     sources   - Cell array of source structures (e.g. {src1, src2, ...}).
%                 Each source struct must contain a .fraction field (scalars,
%                 summing to 1 across all sources).
%     Ntotal    - Total number of particles to generate across all sources.
%     d         - Problem dimension (2 or 3).
%     acoustics - True for acoustic simulation, false for elastic.
%
%   OUTPUT:
%     P - Combined particle structure containing all Ntotal particles, with
%         P.x, P.dir, P.perp, P.p, P.t, P.alive, P.dir0 vertically concatenated,
%         plus P.sourceIdx indicating source attribution.
%
%   PARTICLE CONSERVATION:
%     Ni = round(fraction_i * Ntotal), with the remainder absorbed by the
%     last source so that sum(Ni) == Ntotal exactly.
%
%   See also INITIALIZEPARTICLE, RADIATIVETRANSFER.

if ~iscell(sources)
    error('sources must be a CELL array of source structures.');
end

Ns = numel(sources);
if Ns == 0
    error('sources cell array cannot be empty.');
end

% Extract fractions
fractions = zeros(Ns, 1);
for i = 1:Ns
    if ~isfield(sources{i}, 'fraction') || isempty(sources{i}.fraction)
        error('Source %d is missing the mandatory .fraction field.', i);
    end
    fractions(i) = sources{i}.fraction;
end

% Normalize fractions if they do not sum to 1
sumFrac = sum(fractions);
if abs(sumFrac - 1) > 1e-6
    warning('Source fractions summed to %.4g; normalizing to 1.', sumFrac);
    fractions = fractions / sumFrac;
end

% Compute exact particle count per source
counts = round(fractions * Ntotal);
counts(end) = Ntotal - sum(counts(1:end-1)); % Absorb any rounding remainder

% Pre-allocate cells for concatenation
P_cells = cell(Ns, 1);
src_idx_cells = cell(Ns, 1);

for i = 1:Ns
    ni = counts(i);
    if ni <= 0
        continue;
    end
    src_i = sources{i};
    Pi = initializeParticle(ni, d, acoustics, src_i);
    P_cells{i} = Pi;
    src_idx_cells{i} = repmat(uint32(i), ni, 1);
end

% Filter non-empty sources
valid = ~cellfun(@isempty, P_cells);
P_cells = P_cells(valid);
src_idx_cells = src_idx_cells(valid);

% Combine all fields
P.d = d;
P.N = Ntotal;

f_x     = cellfun(@(p) p.x,     P_cells, 'UniformOutput', false); P.x     = vertcat(f_x{:});
f_dir   = cellfun(@(p) p.dir,   P_cells, 'UniformOutput', false); P.dir   = vertcat(f_dir{:});
f_perp  = cellfun(@(p) p.perp,  P_cells, 'UniformOutput', false); P.perp  = vertcat(f_perp{:});
f_p     = cellfun(@(p) p.p,     P_cells, 'UniformOutput', false); P.p     = vertcat(f_p{:});
f_t     = cellfun(@(p) p.t,     P_cells, 'UniformOutput', false); P.t     = vertcat(f_t{:});
f_alive = cellfun(@(p) p.alive, P_cells, 'UniformOutput', false); P.alive = vertcat(f_alive{:});
f_dir0  = cellfun(@(p) p.dir0,  P_cells, 'UniformOutput', false); P.dir0  = vertcat(f_dir0{:});
P.sourceIdx = vertcat(src_idx_cells{:});

end
