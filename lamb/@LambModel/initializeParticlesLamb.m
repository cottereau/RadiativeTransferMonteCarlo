function P = initializeParticlesLamb(M, src, N, opts)
%INITIALIZEPARTICLESLAMB  K1: emission of a point source near an interface.
%   P = initializeParticlesLamb(M, src, N)
%   P = initializeParticlesLamb(M, src, N, 'xy', [x0 y0], 'ZetaBand', 0.8:0.02:1.2)
%
%   M    model from lambRTEModel;  src  lamb.PointSource (force or moment, depth src.zs,
%        z > 0 in the lower solid, z < 0 in the fluid);  N  number of particles.
%   ZetaBand  optional relative spread of the source depth in wavelengths, i.e. of the
%        frequency (zeta = k_beta1*z scales with omega): the particles are split into equal
%        batches, one per value, each weighted by its own radiated power. This is the
%        band-averaged emission of decision D4.
%
%   |zs| <  M.dsrc : Lamb emission. The channel c is drawn with probability W_c/W_tot from
%                    lamb.rte.channelQ; body directions (k, phi) are drawn jointly from the
%                    node density c_n |d_n R(-phi) s|^2 (k uniform in the mapped variable of
%                    the node's quadrature cell, phi in M.phiBins with jitter). The particles head AWAY from the interface: the
%                    pattern already contains the coherent interface reflection.
%                    Interface-wave particles start on the interface with azimuth drawn from
%                    |h_j R(-phi) s|^2. Leaky particles are drawn from the Pf nodes under the
%                    Lorentzian peak of each near-real leaky pole.
%   |zs| >= M.dsrc : whole-space emission (lower solid only); the particles going toward the
%                    interface later cross it with Zoeppritz (propagateParticleLamb).
%   The path is chosen by the nominal src.zs, so a ZetaBand never mixes the two.
%   Positions: the source point; medium-2 body particles and interface-wave particles start
%   on the interface above the source. Error O(|zs|) <= dsrc, i.e. a few wavelengths (R7).
%
%   P fields (N-by-...): .x .dir .perp .mode (1 P, 2 S, 10+j interface wave j, 20+l leaky l)
%   .side (1: z > 0, 2: z < 0) .t .E (energy per particle, sum = radiated power) .alive
arguments
    M LambModel
    src (1,1) lamb.PointSource
    N (1,1) double {mustBeInteger, mustBePositive}
    opts.xy (1,2) double = [0 0]
    opts.ZetaBand double = 1
end
nb = numel(opts.ZetaBand);
Nb = diff(round(linspace(0, N, nb + 1)));
parts = cell(1, nb);
far = abs(src.zs) >= M.dsrc;                         % path chosen by the NOMINAL depth
for b = 1:nb
    zs = src.zs*opts.ZetaBand(b);
    parts{b} = emit(M, src, zs, Nb(b), opts.xy, far);
    parts{b}.E = parts{b}.E/nb;                      % band average of the radiated power
end
P = parts{1};
for b = 2:nb
    for f = fieldnames(P).'
        P.(f{1}) = [P.(f{1}); parts{b}.(f{1})];
    end
end
P.N = size(P.x, 1);
end

% ------------------------------------------------------------------------------------------
function P = emit(M, src, zs, N, xy, far)
s = lamb.rte.sourceVector(src);
zeta = M.kb1*zs;
if strcmp(src.kind, 'moment') && abs(zeta) < M.zetaMin
    zeta = sign(zeta + (zeta == 0))*M.zetaMin;       % R11: no point moment ON the interface
end
if far
    if zs < 0, error('initializeParticlesLamb:far', 'Far sources must be in the lower solid.'); end
    T = M.Tws;  zeta = 1;                            % whole space: depth only sets phases
else
    T = M.T;
end
[~, W] = lamb.rte.channelQ(T, zeta, M.omega, s);

% channel list with powers: body (minus leaky peaks), interface waves, leaky
names = fieldnames(W.body).';
pw = cellfun(@(n) W.body.(n), names);
iPf = find(strcmp(names, 'Pf'));
if ~isempty(iPf), pw(iPf) = pw(iPf) - sum(W.leaky); end
nIW = numel(W.iw);  nL = numel(W.leaky);
pw = [pw, W.iw, W.leaky];
kind = [repmat({'body'}, 1, numel(names)), repmat({'iw'}, 1, nIW), repmat({'leaky'}, 1, nL)];
idx  = [1:numel(names), 1:nIW, 1:nL];
pw = max(pw, 0);
Wtot = W.total;
cnt = accumarray(sampleDiscrete(pw(:), N), 1, [numel(pw) 1]).';   % multinomial counts

P = blank(0);
for c = find(cnt > 0)
    n = cnt(c);
    switch kind{c}
        case 'body'
            name = names{idx(c)};
            mask = [];
            if strcmp(name, 'Pf') && nL > 0
                mask = ~any(channelMask(T, zeta, name), 2);   % exclude the leaky peaks
            end
            Q = sampleBody(M, T, name, zeta, s, n, mask, far, zs, xy);
        case 'iw'
            Q = sampleIW(M, T, idx(c), zeta, s, n, xy);
        case 'leaky'
            mask = channelMask(T, zeta, 'Pf');
            Q = sampleBody(M, T, 'Pf', zeta, s, n, mask(:, idx(c)), false, zs, xy);
            Q.mode(:) = 20 + idx(c);
            Q.side(:) = 1;  Q.x(:, 3) = 0;
            Q.dir = [Q.dir(:, 1:2)./vecnorm(Q.dir(:, 1:2), 2, 2), zeros(n, 1)];
    end
    P = cat1(P, Q);
end
P.E = Wtot/N*ones(size(P.x, 1), 1);
end

function Q = sampleBody(M, T, name, zeta, s, n, mask, far, zs, xy)
[d, c, k] = lamb.rte.channelAmplitudes(T, name, zeta, M.omega);
if ~isempty(mask), c = c.*mask; end
A = d*reshape(pagemtimes(M.Rneg, s), 9, []);            % N-by-NPhiBins amplitudes
w = c.*abs(A).^2;
[in, im] = ind2sub(size(w), sampleDiscrete(w(:), n));
% k inside the node's own quadrature cell (uniform in the mapped variable), phi jittered
% inside its bin
chn = M.chan(strcmp({M.chan.name}, name));
if far, chn = M.chan(strcmp({M.chan.name}, [name(1:end-1) '1'])); end
km = M.omega/chn.c;
sd = T.solid;  if zeta < 0, sd = T.fluid; end
cells = sd.body(strcmp({sd.body.name}, name)).cells;     % at omega0: scale k to omega
kk  = lamb.quad.cellSample(cells(in, :), rand(n, 1))*(M.omega/T.omega0);
dph = M.phiBins(2) - M.phiBins(1);
ph  = M.phiBins(im).' + (rand(n, 1) - 0.5)*dph;
sinT = min(kk/km, 1);  cosT = sqrt(1 - sinT.^2);
side = chn.side;
down = (side == 1);                                      % away from the interface
if far, down = (name(end) == '1'); end                   % whole space: '2' channels go up
dir = [sinT.*cos(ph), sinT.*sin(ph), (2*down - 1)*cosT];
l = [cos(ph), sin(ph), zeros(n, 1)];  t = [-sin(ph), cos(ph), zeros(n, 1)];
switch name(1:2)
    case 'SV', perp = cross(dir, t, 2);
    case 'SH', perp = t;
    otherwise, perp = randPerp(dir);
end
Q = blank(n);
Q.dir = dir;  Q.perp = perp;  Q.mode(:) = chn.mode;
if far
    Q.side(:) = 1;  Q.x = repmat([xy zs], n, 1);
else
    Q.side(:) = side;
    if side == 1 && zs >= 0 || side == 2 && zs < 0
        Q.x = repmat([xy zs], n, 1);                     % the source's own medium
    else
        Q.x = repmat([xy (3 - 2*side)*1e-9*abs(M.dsrc)], n, 1);   % just across (side 1: z > 0)
    end
end
end

function Q = sampleIW(M, T, j, zeta, s, n, xy)
h = lamb.rte.channelAmplitudes(T, sprintf('IW%d', j), zeta, M.omega);
w = abs(h*reshape(pagemtimes(M.Rneg, s), 9, [])).^2;
im = sampleDiscrete(w(:), n);
dph = M.phiBins(2) - M.phiBins(1);
ph = M.phiBins(im).' + (rand(n, 1) - 0.5)*dph;
Q = blank(n);
Q.x = repmat([xy 0], n, 1);
Q.dir = [cos(ph), sin(ph), zeros(n, 1)];
Q.perp = repmat([0 0 1], n, 1);
Q.mode(:) = 10 + j;  Q.side(:) = 1;
end

function i = sampleDiscrete(w, n)
% n independent draws of an index with probability w(i)/sum(w) (w >= 0)
nz  = find(w > 0);
cdf = cumsum(w(nz));
ok  = [true; diff(cdf) > 0];                    % rounding can repeat an edge: merge
nz  = nz(ok);  cdf = cdf(ok);
[~, ~, b] = histcounts(rand(n, 1)*cdf(end), [0; cdf]);
i = nz(max(b, 1));
end

function m = channelMask(T, zeta, name)
sd = T.solid;  if zeta < 0, sd = T.fluid; end
m = sd.body(strcmp({sd.body.name}, name)).leakMask;
end

function P = blank(n)
P = struct('x', zeros(n, 3), 'dir', zeros(n, 3), 'perp', zeros(n, 3), 'mode', zeros(n, 1), ...
           'side', zeros(n, 1), 't', zeros(n, 1), 'E', zeros(n, 1), 'alive', true(n, 1));
end

function P = cat1(P, Q)
for f = fieldnames(P).'
    P.(f{1}) = [P.(f{1}); Q.(f{1})];
end
end

function p = randPerp(d)
p = cross(d, randn(size(d)), 2);
p = p./vecnorm(p, 2, 2);
end
