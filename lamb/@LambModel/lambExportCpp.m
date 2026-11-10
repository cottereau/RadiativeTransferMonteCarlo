function lambExportCpp(M, file, P)
%LAMBEXPORTCPP  Write a free-surface Lamb model (and optionally particles) to HDF5 for the C++ engine.
%   lambExportCpp(M, file)        M from lambRTEModel (+ lambScatteringModel, + .geom, .edge)
%   lambExportCpp(M, file, P)     also the particles P (initializeParticlesLamb), group /P
%
%   Layout: MATLAB column-major
%   arrays as they are (C++ indexes them Fortran-style, first index fastest); complex arrays as
%   <name>_re and <name>_im; logical as uint8; scalars as 1x1. Only what the transport needs at
%   run time is written. A random medium (het.Phi) is written as the rational fit of lambDepthCorr
%   (Shat(q) = sum_l r_l/(q^2 - s_l)), which is also what g(K, y, dk) uses.
arguments
    M LambModel
    file char
    P struct = struct()
end
I = M.iface;
assert(isempty(I.upper), 'lambExportCpp: free surface only');
if isfile(file), delete(file); end
L = I.lower;
put(file, '/medium/alpha', L.alpha);  put(file, '/medium/beta', L.beta);  put(file, '/medium/rho', L.rho);
put(file, '/medium/lambda', L.lambda);  put(file, '/medium/mu', L.mu);
put(file, '/omega', M.omega);  put(file, '/kb1', M.kb1);  put(file, '/dsrc', M.dsrc);
put(file, '/dslab', M.dslab);  put(file, '/zetaMin', M.zetaMin);
put(file, '/iw/c', [M.iw.c]);
% --- half-space tables (solid side) ------------------------------------------------------------
T = M.T;  put(file, '/T/omega0', T.omega0);  put(file, '/T/beta', T.iface.lower.beta);
if ~isempty(M.zMirror), put(file, '/zMirror', M.zMirror); end
for c = 1:numel(T.solid.body)
    B = T.solid.body(c);  g = sprintf('/T/body%d/', c);
    put(file, [g 'type'], find(strcmp(B.name, {'P1', 'SV1', 'SH1'})));
    for f = {'k', 'cells', 'c', 'Tr', 'tow', 'nuTow', 'away', 'nuAway'}, put(file, [g f{1}], B.(f{1})); end
end
put(file, '/T/nbody', numel(T.solid.body));
for j = 1:numel(T.solid.iw)
    Pj = T.solid.iw(j);  g = sprintf('/T/iw%d/', j);
    for f = {'k', 'fac', 'q', 'tow', 'nuTow'}, put(file, [g f{1}], Pj.(f{1})); end
end
put(file, '/T/niw', numel(T.solid.iw));
% --- scattering ------------------------------------------------------------------------------
if ~isempty(M.het)
    het = M.het{1};  psd = isfield(het, 'law');
    put(file, '/het/psd', double(psd));
    if psd
        [~, R] = lambDepthCorr(het, 0, 1, 0);
        put(file, '/het/a', het.a);  put(file, '/het/cv', het.cv);  put(file, '/het/corr', het.corr);
        put(file, '/het/s', R.s);  put(file, '/het/r', R.r);
    else
        for f = {'n', 'V', 'drho', 'dlam', 'dmu'}, put(file, ['/het/' f{1}], het.(f{1})); end
    end
    put(file, '/bulk/sigma', M.bulk{1}.sigma);  put(file, '/bulk/pP', M.bulk{1}.pP);  put(file, '/bulk/rate', M.bulk{1}.rate);
    for j = 1:numel(M.iwT)
        X = M.iwT(j);  g = sprintf('/iwT%d/', j);
        for f = {'z', 'wz', 'nsig', 'mu'}, put(file, [g f{1}], X.(f{1})); end
        if ~isempty(X.psd)                                          % stacked over the depth nodes
            put(file, [g 'psd/Hc'], cat(3, X.psd.Hc));  put(file, [g 'psd/s'], cat(4, X.psd.s));
            put(file, [g 'psd/K'], cat(1, X.psd.K));    put(file, [g 'psd/H9'], cat(3, X.psd.H9));
        end
    end
    if ~isempty(M.loc)
        Lc = M.loc;  g = '/loc/';
        put(file, [g 'ctype'], cellfun(@(t) find(strcmp(regexprep(t, '[12f]$', ''), {'P', 'SV', 'SH'})), Lc.types));
        for f = {'c', 'k', 'z', 'dz', 'Hin', 'Hout', 'rin', 'rout', 'f', 'fin', 'fout', 'CinA', 'CoutA'}
            put(file, [g f{1}], Lc.(f{1}));
        end
        if isfield(Lc, 'fnum'), put(file, [g 'fnum'], Lc.fnum);  put(file, [g 'fden'], Lc.fden); end
        if isfield(Lc, 'psd') && ~isempty(Lc.psd)
            X = Lc.psd;
            for f = {'Hcin', 'Hcout', 'G', 'sref', 'K', 'zref', 'c', 'sm', 'w'}, put(file, [g 'psd/' f{1}], X.(f{1})); end
        end
    end
end
% --- geometry and edges --------------------------------------------------------------------------
if ~isempty(M.geom), put(file, '/geom/R', M.geom.R);  put(file, '/geom/h', M.geom.h); end
if ~isempty(M.edge)
    E = M.edge;
    for f = {'psi', 'T', 'R', 'P', 'S', 'thEdges', 'pP', 'pS'}, put(file, ['/edge/' f{1}], E.(f{1})); end
end
% --- particles -----------------------------------------------------------------------------------
if isfield(P, 'x')
    for f = {'x', 'dir', 'perp', 'mode', 'side', 't', 'E'}, put(file, ['/P/' f{1}], P.(f{1})); end
    if isfield(P, 'face'), put(file, '/P/face', P.face); end
end
end

function put(file, name, v)
if islogical(v), v = uint8(v); end
if ~isreal(v)
    put(file, [name '_re'], real(v));  put(file, [name '_im'], imag(v));  return
end
if isempty(v), v = zeros(1, 0); end
sz = size(v);  if isvector(v) && numel(sz) == 2, sz = [numel(v) 1]; v = v(:); end
if any(sz == 0), sz(sz == 0) = 1;  v = zeros(sz); name = [name '_empty']; end
h5create(file, name, sz, 'Datatype', class(v));
h5write(file, name, v);
end
