function F = lambRandomForms(Gr, M, rule)
%LAMBRANDOMFORMS  Rule-B forms over the components of lambRandomGram.
%   F = lambRandomForms(Gr, M)
%
%   lambLocalForms in the component basis of Gr (vectors x over the components; a point source
%   S x): every 9-by-9 form S H S' of lambLocalForms is H here (nu-by-nu), and its removal value
%   tr(Q_t S H S') becomes real(sum(H .* Gr.G, 'all')) (random medium: the nonlocal Gram).
%     F.Hin{b}, F.Hout{b}  in / out forms of the propagating types b of the medium of Gr.z
%     rule (3rd argument): 'flip' (default, rule B': f = 1/2, flips F.fin in -> out, F.fout
%                          out -> in, see lambLocalForms) or 'B'
%     F.f(b)               evanescent share (R_in = R_out); f = F.fnum/F.fden, both smooth in k
%                          (f has poles where the evanescent own power crosses zero: interpolate
%                          fnum and fden, not f)
%     F.rin(b), F.rout(b)  removal values (per d2k and per F n, as lambStationarity)
%     F.Htun{a}, F.rtun(a) tunnelling forms / values of the OTHER medium's in-cells
%   The emission pattern of a form H is sum_{uv} H_uv T_u conj(T_v) w_uv (lambRandomGram's pair
%   weights), split by the coherent labels.
if nargin < 3, rule = 'flip'; end                                  % rule B' (see lambLocalForms)
val = @(H) real(sum(H.*Gr.G, 'all'));
sym = @(x, y) (x*y' + y*x')/2;
nq = numel(Gr.types);  nu = Gr.nu;
so = cell(1, nq);  for q = 1:nq, so{q} = zeros(nu, 1);  so{q}(Gr.iOut(q)) = 1; end
st = Gr.states;  m = Gr.m;  mine = find([st.med] == m);
F.f = nan(1, nq);  F.rin = zeros(1, nq);  F.rout = F.rin;  F.Hin = cell(1, nq);  F.Hout = F.Hin;
F.fnum = nan(1, nq);  F.fden = nan(1, nq);  F.Ain = zeros(1, nq);  F.Aout = F.Ain;  F.O = F.Ain;  F.fin = F.Ain;  F.fout = F.Ain;
for b = find(Gr.prop)
    i = mine([st(mine).type] == b);
    if isempty(i), continue, end                                   % grazing within 1e-12
    si = st(i).si;  se = st(i).se;  c = st(i).c;  Nb = st(i).N;
    Hown = se*se' + 2*sym(si, se);  Hin = si*si';
    for q = find(Gr.prop)
        Hown = Hown + sym(conj(c(q))*se, so{q});  Hin = Hin + sym(conj(c(q))*si, so{q});
    end
    Hout = so{b}*so{b}'/Nb;
    for j = 1:numel(st)
        Hout = Hout + sym(conj(st(j).c(b))*(st(j).si + st(j).se), so{b})/st(j).N;
    end
    ain = val(Hin);  aout = val(Hout);  ob = val(Hown);
    if strcmp(rule, 'flip'), f = 0.5;
    elseif abs(ob) > 1e-12*max(abs([ain, aout*Nb, realmin])), f = (aout*Nb - ain + ob)/(2*ob); else, f = 0.5;
    end
    F.Hin{b} = (Hin + f*Hown)/Nb;  F.Hout{b} = Hout + (1 - f)*Hown/Nb;  F.f(b) = f;
    F.fnum(b) = aout*Nb - ain + ob;  F.fden(b) = 2*ob;
    F.Ain(b) = ain/Nb;  F.Aout(b) = aout;  F.O(b) = ob/Nb;
    if strcmp(rule, 'flip'), D = ain/Nb - aout;  F.fin(b) = max(-D, 0);  F.fout(b) = max(D, 0);  F.fnum(b) = 0.5;  F.fden(b) = 1; end
    F.rin(b) = val(F.Hin{b});  F.rout(b) = val(F.Hout{b});
end
cells = lambMediumCells(M.iface);  F.Htun = {};  F.rtun = [];
if numel(cells) == 2
    no = numel(cells(3 - m).types);  F.Htun = cell(1, no);  F.rtun = zeros(1, no);
    for j = find([st.med] ~= m)
        H = st(j).se*st(j).se';
        for q = find(Gr.prop), H = H + sym(conj(st(j).c(q))*st(j).se, so{q}); end
        F.Htun{st(j).type} = H/st(j).N;  F.rtun(st(j).type) = val(F.Htun{st(j).type});
    end
end
end
