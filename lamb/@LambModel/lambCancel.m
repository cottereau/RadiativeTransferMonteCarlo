function P = lambCancel(M, P, opts)
%LAMBCANCEL  Weight cancellation and population control for signed particles.
%   P = lambCancel(M, P)
%   P = lambCancel(M, P, 'NTheta', 12, 'DzSlab', d/6, 'DzBulk', 2*pi, 'Target', N0)
%
%   The local slab transport emits signed particles (the coherent near-interface field is not
%   a positive mixture of independent plane waves); without control their |E| grows by 8-20 %
%   per scattering event. Once per step, alive particles are grouped into phase-space cells:
%     body : medium side, wave class (P | S mostly SV | S mostly SH, by aSV^2 > 1/2), sense
%            (toward / away from z = 0), theta bin (NTheta, uniform in theta), depth bin
%            (DzSlab inside |z| < d, DzBulk beyond), optional azimuth (NPhi) and horizontal
%            position (Dxy) bins;
%     IW / leaky : mode, optional azimuth and horizontal bins.
%   In every cell with both signs the particles are replaced by n = max(1, round(count
%   |sum E| / sum |E|)) particles drawn (with replacement, ~|E|) from the majority sign, each
%   carrying (sum E)/n: the signed energy of every cell is kept exactly, the distribution inside
%   the cell is represented by its majority-sign particles (bias at the cell scale). Then the
%   population is brought back to Target (default: the count at the first call): down by Russian
%   roulette when above 1.2 Target (|E| < w0 survives with probability |E|/w0 at energy +-w0),
%   up by splitting the largest-|E| particles in two (dead rows are removed).
%   Defaults resolve what gate G4.8 balances (12 theta bins, lambda/2 windows) and ignore
%   azimuth and horizontal position (horizontally invariant, isotropic problems such as G4.4).
arguments
    M LambModel
    P struct
    opts.NTheta (1,1) double = 12
    opts.DzSlab (1,1) double = NaN
    opts.DzBulk (1,1) double = 2*pi
    opts.NPhi   (1,1) double = 1
    opts.Dxy    (1,1) double = Inf
    opts.Target (1,1) double = NaN
end
if isnan(opts.DzSlab), opts.DzSlab = M.dslab/6; end
if ~isfield(P, 'N0'), P.N0 = P.N; end
if isnan(opts.Target), opts.Target = P.N0; end
a = find(P.alive(:));  n = numel(a);
% --- cell keys --------------------------------------------------------------------------------
body = P.mode(a) <= 2;
cls = P.mode(a);                                                   % 1 P, 2 S, >= 10 IW / leaky
s = find(body & P.mode(a) == 2);
if ~isempty(s)
    L = lambConventionMap(P.dir(a(s), :), P.perp(a(s), :), [0 0 1]);
    cls(s) = 2 + (L.aSV.^2 < 0.5);                                 % 2 SV-like, 3 SH-like
end
z = P.x(a, 3);  az = abs(z);
toward = (P.side(a) == 1 & P.dir(a, 3) < 0) | (P.side(a) == 2 & P.dir(a, 3) > 0);
th = acos(min(abs(P.dir(a, 3)), 1));
tb = min(floor(th/(pi/2)*opts.NTheta), opts.NTheta - 1);
d = M.dslab;
zb = floor(az/opts.DzSlab);
out = az >= d;  zb(out) = ceil(d/opts.DzSlab) + floor((az(out) - d)/opts.DzBulk);
ph = atan2(P.dir(a, 2), P.dir(a, 1));
pb = min(floor(mod(ph, 2*pi)/(2*pi)*opts.NPhi), opts.NPhi - 1);
if isfinite(opts.Dxy), xb = floor(P.x(a, 1:2)/opts.Dxy); else, xb = zeros(n, 2); end
key = [P.side(a), cls, toward, tb, zb, pb, xb];
key(~body, [1 3 4 5]) = 0;                                         % IW / leaky: mode (+ phi, xy)
[~, ~, cell] = unique(key, 'rows');
% --- cancellation -------------------------------------------------------------------------------
E = P.E(a);
Sn = accumarray(cell, E);  Sa = accumarray(cell, abs(E));  cnt = accumarray(cell, 1);
mixed = find(accumarray(cell, E < 0) > 0 & accumarray(cell, E > 0) > 0);
keepRows = true(n, 1);  newSrc = zeros(0, 1);  newE = zeros(0, 1);
for c = mixed.'
    in = find(cell == c);  net = Sn(c);
    keepRows(in) = false;
    if net == 0, continue, end
    maj = in(sign(E(in)) == sign(net));
    nn = max(1, round(cnt(c)*abs(net)/Sa(c)));
    cw = cumsum(abs(E(maj)));
    pick = maj(sum(cw < rand(nn, 1).'*cw(end), 1) + 1);
    newSrc = [newSrc; pick(:)]; %#ok<AGROW>
    newE = [newE; repmat(net/nn, nn, 1)]; %#ok<AGROW>
end
rows = [a(keepRows); a(newSrc)];
Ekeep = [E(keepRows); newE];
P = takeRows(P, rows);
P.E = Ekeep;
% --- population control: split the largest |E| back to the target count ---------------------
short = round(opts.Target) - P.N;
if P.N > 1.2*opts.Target                                           % Russian roulette down to Target
    a = abs(P.E);  lo = 0;  hi = max(a);
    for it = 1:60                                                  % w0: sum(min(1, |E|/w0)) = Target
        w0 = (lo + hi)/2;
        if sum(min(1, a/w0)) > opts.Target, lo = w0; else, hi = w0; end
    end
    w0 = hi;  keep = rand(P.N, 1) < min(1, a/w0);
    sm = keep & a < w0;  P.E(sm) = sign(P.E(sm))*w0;               % survivors keep E in expectation
    P = takeRows(P, find(keep));
end
short = round(opts.Target) - P.N;
while short > 0 && P.N > 0                                         % rounds: each row split once
    [~, o] = sort(abs(P.E), 'descend');
    sp = o(1:min(short, P.N));
    P.E(sp) = P.E(sp)/2;
    P = takeRows(P, [(1:P.N).'; sp(:)]);
    short = round(opts.Target) - P.N;
end
end

function P = takeRows(P, rows)
% keep / duplicate the given rows of every per-particle field
n0 = P.N;
for f = fieldnames(P).'
    F = P.(f{1});
    if strcmp(f{1}, 'N') || strcmp(f{1}, 'N0') || size(F, 1) ~= n0, continue, end
    P.(f{1}) = F(rows, :);
end
P.N = numel(rows);
P.alive = true(P.N, 1);
end
