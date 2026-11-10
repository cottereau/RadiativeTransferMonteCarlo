function loc = lambLocalTables(M, opts)
%LAMBLOCALTABLES  Tables of the local slab transport (rule B).
%   loc = lambLocalTables(M)
%   loc = lambLocalTables(M, 'Preset', 'fine')            % 'coarse' | 'standard' (default) | 'fine'
%   loc = lambLocalTables(M, 'NkSeg', 40, 'Nz', 128)      % explicit resolution (overrides the preset)
%   Resolution presets (NkSeg / Nz), measured off the nodes against lambLocalForms
%   (test_lamb_G61_tables, 2026-10-06; criteria: local rate within 1 % of its peak, layer-integrated
%   scattered power per channel within 0.5 %):
%     'coarse'    40 / 128   rate 2.5-4.8 %, layer 1.6-3.8 % (worst near grazing / critical angles)
%     'standard'  80 / 256   rate 0.9-2.3 %, layer 1.1-1.7 % (0.1-0.3 % away from grazing/critical)
%     'fine'     160 / 512   rate <= 0.47 %, layer <= 0.33 % (free surface, solid-fluid dmu, Al/W dmu)
%   Table time grows fast: Al/W dmu ~70 s / 16 min / 1 h 51 for the three presets.
%
%   M needs lambScatteringModel's .het and .dslab. For each medium m (1 lower, 2 upper) on a
%   grid k (common to both media, breakpoints at every w/c, sine-clustered in each segment)
%   x z (slab depths of m, |z| = d (i - 1/2)/Nz, about 10 nodes per interference fringe):
%     loc(m).types, .c, .k (1-by-Nk), .z (1-by-Nz, signed), .dz
%     .Hin, .Hout   9-by-9-by-Nk-by-Nz-by-ntypes forms of lambLocalForms (emission patterns;
%                   single precision: only the shape of the pattern is taken from them)
%     .rin, .rout   Nk-by-Nz-by-ntypes local rate per unit path: n trace(Qt H) cos(theta_b)
%                   (0 where type b does not propagate)
%     .Cin, .Cout   Nk-by-(Nz+1)-by-ntypes cumulative int_0^|z| rate d|z| at |z| = (0:Nz) dz
%     .CinA, .CoutA the same for |rate| (event sampling; local rates can be negative)
%     .f            Nk-by-Nz-by-ntypes evanescent share (emission of same-medium channels; 1/2 for
%                   rule B', the default)
%     .fin, .fout   Nk-by-Nz-by-ntypes rule-B' mirror-flip rates per unit path (in -> out, out -> in);
%                   CinA / CoutA include them (event = scattering or flip)
%     .Htun         9-by-9-by-Nk-by-Nz-by-ntypesO tunnelling forms owned by the in-cells of the
%                   OTHER medium, events at depths of m
%     .qtun         Nk-by-Nz-by-ntypesO  n trace(Qt Htun) per unit depth; .Ptun Nk-by-ntypesO
%                   the per-arrival probability density integral int qtun d|z| (use
%                   1 - exp(-Ptun)); .Ctun cumulative, for sampling the event depth;
%                   .PtunA, .CtunA the same for |qtun|
%   Random medium (M.het{m}.law set): rates, f and qtun from the nonlocal rule B
%   (lambRandomGram / lambRandomForms, n = 1); Hin / Hout / Htun are then PROPOSAL forms (the
%   point limit S Hc S' of the same component forms) and L.psd holds what the emission needs to
%   reweight to the random medium (lambEmitLocal): component slots 1..nt incident of type b,
%   nt+1..2nt outgoing of type q;  .Hcin .Hcout (2nt x 2nt x Nk x Nz x nt), .Hctun, .G (Gram,
%   2nt x 2nt x Nk x Nz), .sref (9 x 2nt x 3 x Nk, sources at z = .zref), .K (2nt x Nk, s ~
%   e^{i K z}), .het.
%   Checks printed by the caller (see testLambPhase4 G4.8/G4.9): interpolation vs direct forms.
arguments
    M LambModel
    opts.Preset char {mustBeMember(opts.Preset, {'coarse', 'standard', 'fine'})} = 'standard'
    opts.NkSeg (1,1) double = NaN            % NaN: from the preset
    opts.Nz    (1,1) double = NaN
    opts.Rule  char {mustBeMember(opts.Rule, {'flip', 'B'})} = 'flip'   % rule B' (default) or B
end
pre = struct('coarse', [40 128], 'standard', [80 256], 'fine', [160 512]);
if isnan(opts.NkSeg), opts.NkSeg = pre.(opts.Preset)(1); end
if isnan(opts.Nz), opts.Nz = pre.(opts.Preset)(2); end
I = M.iface;  w = M.omega;  d = M.dslab;
cells = lambMediumCells(I);  nm = numel(cells);
allc = [cells.c];
kb = unique(w./allc);  kb = [0, kb(:).'];
t = ((1:opts.NkSeg) - 0.5)/opts.NkSeg;
k = [];
for s = 1:numel(kb) - 1
    k = [k, kb(s) + (kb(s+1) - kb(s))*(1 + sin(pi*(t - 0.5)))/2]; %#ok<AGROW>
end
Nk = numel(k);
for m = 1:nm
    sg = 1 - 2*(m == 2);
    zq = sg*d*((1:opts.Nz) - 0.5)/opts.Nz;
    nt = numel(cells(m).types);  nto = 0;  if nm == 2, nto = numel(cells(3 - m).types); end
    L = struct('types', {cells(m).types}, 'c', cells(m).c, 'k', k, 'z', zq, 'dz', d/opts.Nz);
    L.Hin = complex(zeros(9, 9, Nk, opts.Nz, nt, 'single'));  L.Hout = L.Hin;
    L.rin = zeros(Nk, opts.Nz, nt);  L.rout = L.rin;  L.f = 0.5*ones(Nk, opts.Nz, nt);
    L.fin = zeros(Nk, opts.Nz, nt);  L.fout = L.fin;               % rule B' flip rates per unit path
    L.Htun = complex(zeros(9, 9, Nk, opts.Nz, nto, 'single'));  L.qtun = zeros(Nk, opts.Nz, nto);
    if isfield(M.het{m}, 'law')
        L = randomTables(M, m, L, k, zq, nt, nto, opts.Rule);
    else
    n0 = M.het{m}.n;
    for iz = 1:opts.Nz
        Qt = lambPowerForm(M, zq(iz));
        for ik = 1:Nk
            F = lambLocalForms(M, m, k(ik), zq(iz), Qt, opts.Rule);
            for b = find(F.prop)
                L.Hin(:, :, ik, iz, b) = F.Hin{b};  L.Hout(:, :, ik, iz, b) = F.Hout{b};
                L.rin(ik, iz, b) = n0*real(trace(Qt*F.Hin{b}))*F.cos(b);
                L.rout(ik, iz, b) = n0*real(trace(Qt*F.Hout{b}))*F.cos(b);
                L.f(ik, iz, b) = F.f(b);
                L.fin(ik, iz, b) = n0*F.fin(b)*F.cos(b);  L.fout(ik, iz, b) = n0*F.fout(b)*F.cos(b);
            end
            for a = find(F.propO)
                L.Htun(:, :, ik, iz, a) = F.Htun{a};
                L.qtun(ik, iz, a) = n0*real(trace(Qt*F.Htun{a}));
            end
        end
    end
    end
    cum = @(r) cat(2, zeros(Nk, 1, size(r, 3)), cumsum(r, 2)*L.dz);   % midpoint rule, node 0 at z = 0
    L.Cin = cum(L.rin);  L.Cout = cum(L.rout);  L.Ctun = cum(L.qtun);
    L.CinA = cum(abs(L.rin) + L.fin);  L.CoutA = cum(abs(L.rout) + L.fout);  L.CtunA = cum(abs(L.qtun));   % sampling (events + flips)
    L.Ptun = reshape(L.Ctun(:, end, :), Nk, []);  L.PtunA = reshape(L.CtunA(:, end, :), Nk, []);
    loc(m) = L; %#ok<AGROW>
end
end

function L = randomTables(M, m, L, k, zq, nt, nto, rule)
% random-medium rule B on the (k, z) grid (see the header); ponytail: one Gram per node,
% ~0.1 s each (Nk x Nz ~ 2e4: ~40 min per medium), cache the model if reused
het = M.het{m};  Nk = numel(k);  Nz = numel(zq);  nf = 2*nt;
L.fnum = nan(Nk, Nz, nt);  L.fden = nan(Nk, Nz, nt);             % f = fnum/fden (poles in f)
P = struct('Hcin', complex(zeros(nf, nf, Nk, Nz, nt, 'single')), 'G', complex(zeros(nf, nf, Nk, Nz, 'single')), ...
           'sref', zeros(9, nf, 3, Nk), 'K', zeros(nf, Nk), 'zref', zq(1), 'het', het, ...
           'c', L.c, 'sm', 3 - 2*m, 'w', M.omega);                 % exact K at any k (psdEvent)
P.Hcout = P.Hcin;  P.Hctun = complex(zeros(nf, nf, Nk, Nz, nto, 'single'));
w = M.omega;
for iz = 1:Nz
    for ik = 1:Nk
        Gr = lambRandomGram(M, het, zq(iz), k(ik));
        F = lambRandomForms(Gr, M, rule);
        sl = zeros(1, Gr.nu);                                      % Gr component -> slot
        for b = find(Gr.prop), sl(Gr.iInc(b)) = b; end
        sl(Gr.iOut) = nt + (1:nt);
        pad = @(H) accumarray([repmat(sl(:), Gr.nu, 1), kron(sl(:), ones(Gr.nu, 1))], H(:), [nf nf]);
        H9 = @(H) prop9(Gr.s, H);
        P.G(:, :, ik, iz) = pad(Gr.G);
        if iz == 1, P.sref(:, sl, :, ik) = Gr.s;  P.K(sl, ik) = Gr.K; end
        cs = sqrt(max(1 - (k(ik)*L.c/w).^2, 0));
        for b = find(Gr.prop)
            if isempty(F.Hin{b}), continue, end
            L.Hin(:, :, ik, iz, b) = H9(F.Hin{b});  L.Hout(:, :, ik, iz, b) = H9(F.Hout{b});
            P.Hcin(:, :, ik, iz, b) = pad(F.Hin{b});  P.Hcout(:, :, ik, iz, b) = pad(F.Hout{b});
            L.rin(ik, iz, b) = F.rin(b)*cs(b);  L.rout(ik, iz, b) = F.rout(b)*cs(b);  L.f(ik, iz, b) = F.f(b);
            L.fnum(ik, iz, b) = F.fnum(b);  L.fden(ik, iz, b) = F.fden(b);
            L.fin(ik, iz, b) = F.fin(b)*cs(b);  L.fout(ik, iz, b) = F.fout(b)*cs(b);
        end
        for a = find(~cellfun(@isempty, F.Htun))
            L.Htun(:, :, ik, iz, a) = H9(F.Htun{a});  P.Hctun(:, :, ik, iz, a) = pad(F.Htun{a});
            L.qtun(ik, iz, a) = F.rtun(a);
        end
    end
end
L.psd = P;
end

function H = prop9(s, Hc)
% proposal form sum_j S_j Hc S_j' (point limit of the component form; covariance columns j)
H = zeros(9);
for j = 1:size(s, 3), H = H + s(:, :, j)*Hc*s(:, :, j)'; end
end
