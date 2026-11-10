function F = lambLocalForms(M, m, k, z, Qt, rule)
%LAMBLOCALFORMS  Rule-B scattering forms of the local slab transport at (k, z).
%   F = lambLocalForms(M, m, k, z)
%   F = lambLocalForms(M, m, k, z, Qt)        Qt: total power form at z (else computed)
%
%   m   medium of the depth z (1 lower, z > 0; 2 upper, z < 0);  k  horizontal wavenumber.
%   The diffuse field at z is the incoherent sum of every incident state of both media; each
%   state at z = incident wave s_i (states of m) + evanescent waves s_e + propagating outgoing
%   waves c_q s_q (unit plane waves of m, scatteringState 'Away'); flux normalisation
%   N_a = I_a cos(theta_a). For every type b of m propagating at k (9-by-9 Hermitian forms,
%   physical frame, azimuth 0; per d2k and per F n, i.e. trace(Qt H) is the removal value of
%   lambStationarity; the rate per unit path and scatterer density is trace(Qt H) cos(theta_b)):
%     F.Hin{b}   in-cell  : [s_i s_i' + 1/2 sum_q (c_q' s_i s_q' + c_q s_q s_i') + f Hown] / N_b
%     F.Hout{b}  out-cell : s_b s_b'/N_b + (1 - f) Hown/N_b
%                           + sum_a 1/2 (c_ba' (s_i + s_e)_a s_b' + c_ba s_b (s_i + s_e)_a') / N_a
%     Hown = s_e s_e' + s_i s_e' + s_e s_i' + 1/2 sum_q (c_q' s_e s_q' + c_q s_q s_e')
%     F.f(b)     evanescent share, chosen so that trace(Qt Hin) = trace(Qt Hout)
%   Rule (6th argument, default 'flip' = rule B'): f = 1/2 and the in/out imbalance
%   D = A_in - A_out (A_in = trace(Qt Hin_raw)/N_b, A_out = trace(Qt Hout_raw)) is moved by mirror
%   flips: F.fin(b) = max(-D, 0) (in -> out), F.fout(b) = max(D, 0) (out -> in), per d2k and n as
%   the removal values; removal = trace(Qt Hin) + fin = trace(Qt Hout) + fout. 'B' = the original
%   rule B (f from R_in = R_out; has poles where the evanescent own power crosses zero).
%   F.Ain, F.Aout, F.O (= own evanescent / N_b) for either rule.
%   and for every type a of the OTHER medium propagating at k (tunnelling, owned by its in-cell):
%     F.Htun{a}  [s_e s_e' + 1/2 sum_q (c_q' s_e s_q' + c_q s_q s_e')] / N_a
%   The emission pattern of a form H in the direction with 9-vector amplitude row a is a H a'.
%   F.prop(b), F.cos(b) (medium m), F.propO(a), F.cosO(a) (other medium), F.types, F.typesO.
%
%   See also LAMBSTATIONARITY (gate G4.8, which uses these forms), LAMBLOCALTABLES.
if nargin < 5 || isempty(Qt), Qt = lambPowerForm(M, z); end
if nargin < 6, rule = 'flip'; end
I = M.iface;  w = M.omega;
cells = lambMediumCells(I);
S = allStates(M, m, k, z, cells);
Cm = cells(m);  nq = numel(Cm.types);
F.types = Cm.types;  F.prop = S.prop;  F.cos = sqrt(max(1 - (k*Cm.c/w).^2, 0));
F.Hin = cell(1, nq);  F.Hout = cell(1, nq);  F.f = nan(1, nq);
F.Ain = zeros(1, nq);  F.Aout = F.Ain;  F.O = F.Ain;  F.fin = F.Ain;  F.fout = F.Ain;
sym = @(a, b) (a*b' + b*a')/2;                                     % Hermitian part of a b'
tr = @(H) real(trace(Qt*H));
for b = find(S.prop)
    mi = find(S.med == m & S.type == b);
    Hown = S.se{mi}*S.se{mi}' + 2*sym(S.si{mi}, S.se{mi});
    Hin = S.si{mi}*S.si{mi}';
    for q = find(S.prop)
        Hown = Hown + sym(conj(S.c{mi}(q))*S.se{mi}, S.so{q});
        Hin = Hin + sym(conj(S.c{mi}(q))*S.si{mi}, S.so{q});
    end
    Nb = S.N(mi);
    Hout = S.so{b}*S.so{b}'/Nb;
    for i = 1:numel(S.N)
        Hout = Hout + sym(conj(S.c{i}(b))*(S.si{i} + S.se{i}), S.so{b})/S.N(i);
    end
    ain = tr(Hin);  aout = tr(Hout);  ob = tr(Hown);
    if strcmp(rule, 'flip')                                        % rule B': share 1/2 + flips
        f = 0.5;
    elseif abs(ob) > 1e-12*max(abs([ain, aout*Nb, realmin])), f = (aout*Nb - ain + ob)/(2*ob); else, f = 0.5;
    end
    F.Hin{b} = (Hin + f*Hown)/Nb;
    F.Hout{b} = Hout + (1 - f)*Hown/Nb;
    F.f(b) = f;
    F.Ain(b) = ain/Nb;  F.Aout(b) = aout;  F.O(b) = ob/Nb;
    if strcmp(rule, 'flip')                                        % flips: in -> out, out -> in
        D = ain/Nb - aout;  F.fin(b) = max(-D, 0);  F.fout(b) = max(D, 0);
    end
end
F.typesO = {};  F.Htun = {};  F.propO = [];  F.cosO = [];
if numel(cells) == 2
    Co = cells(3 - m);  F.typesO = Co.types;  F.Htun = cell(1, numel(Co.types));
    F.cosO = sqrt(max(1 - (k*Co.c/w).^2, 0));  F.propO = false(1, numel(Co.types));
    for a = 1:numel(Co.types)
        i = find(S.med == 3 - m & S.type == a);
        if isempty(i), continue, end
        F.propO(a) = true;
        H = S.se{i}*S.se{i}';
        for q = find(S.prop), H = H + sym(conj(S.c{i}(q))*S.se{i}, S.so{q}); end
        F.Htun{a} = H/S.N(i);
    end
end
end

% ------------------------------------------------------------------------------------------
function S = allStates(M, m, k, z, cells)
% every incident state (both media) propagating at k, decomposed at z (medium m) into the
% incident wave s_i (states of m), evanescent waves s_e and outgoing waves c(q) s_o{q};
% the decomposition is done on the fields (u, eps), column-scaled, and checked
I = M.iface;  w = M.omega;  p = k/w;  het = M.het{m};
src = @(St) lamb.rte.equivalentSource(St, het, w);
fld = @(St) [St.u(:); St.eps(:)];
Cm = cells(m);  nq = numel(Cm.types);
S.prop = k < w./Cm.c*(1 - 1e-12);
U = zeros(12, nq);  S.so = cell(1, nq);
for q = 1:nq
    St = lamb.rte.scatteringState(I, Cm.types{q}, p, w, z, 'Mode', 'incident', 'Away', true);
    U(:, q) = fld(St);  S.so{q} = src(St);
end
un = vecnorm(U);
S.si = {};  S.se = {};  S.c = {};  S.N = [];  S.med = [];  S.type = [];
for mm = 1:numel(cells)
    for a = 1:numel(cells(mm).types)
        ca = cells(mm).c(a);
        if k >= w/ca*(1 - 1e-12), continue, end                     % not propagating
        Ss = lamb.rte.scatteringState(I, cells(mm).types{a}, p, w, z, 'Mode', 'state');
        tgt = fld(Ss);  si = zeros(9, 1);
        if mm == m
            Si = lamb.rte.scatteringState(I, cells(mm).types{a}, p, w, z, 'Mode', 'incident');
            tgt = tgt - fld(Si);  si = src(Si);
        end
        c = ((U./un)\tgt)./un.';
        c(abs(c).*un.' < 1e-11*norm(tgt)) = 0;                      % round-off leakage
        assert(norm(U*c - tgt) <= 1e-8*norm(tgt) + 1e-14, 'lambLocalForms: decomposition failed');
        se = zeros(9, 1);
        for q = find(~S.prop), se = se + c(q)*S.so{q}; end
        S.si{end+1} = si;  S.se{end+1} = se;  S.c{end+1} = c;
        S.N(end+1) = Ss.I*sqrt(1 - (k*ca/w)^2);  S.med(end+1) = mm;  S.type(end+1) = a;
    end
end
end
