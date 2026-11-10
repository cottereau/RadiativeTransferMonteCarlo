function [qg,Sg] = structureFactor(arqs,qmin,qmax,nq,ndir)
%% structureFactor - S(q) of a periodic sphere packing
%
%   [q,S] = TMatrix.structureFactor(centerFiles, qmin, qmax, nq, ndir)
%
% S(q) = <|sum_j exp(i q.r_j)|^2>/N, averaged over realizations.
%
% The box is PERIODIC (Lame.py builds the cKDTree with boxsize), so only
% reciprocal-lattice vectors q = 2*pi*(n1/Lx, n2/Ly, n3/Lz) are allowed.
% Sampling random directions evaluates forbidden q and produces a spurious
% peak of order N at low q (previously seen as "S at the peak = 90").
L = [0.20 0.20 0.30];
if ischar(arqs) || isstring(arqs), arqs = cellstr(arqs); end

qg = logspace(log10(qmin),log10(qmax),nq);
edges = [qg(1)*0.9, sqrt(qg(1:end-1).*qg(2:end)), qg(end)*1.1];

nmax = ceil(qmax*max(L)/(2*pi)) + 1;
[i1,i2,i3] = ndgrid(-nmax:nmax,-nmax:nmax,-nmax:nmax);
Q  = [i1(:)/L(1), i2(:)/L(2), i3(:)/L(3)]*2*pi;
qn = vecnorm(Q,2,2);
keep = qn >= edges(1) & qn <= edges(end);
Q = Q(keep,:); qn = qn(keep);

Sg = zeros(1,nq); cont = zeros(1,nq);
rng(7)
for ia = 1:numel(arqs)
    d = readmatrix(arqs{ia},'FileType','text','CommentStyle','#')*1e-4;  % 1 un = 0.1 mm
    d = d(all(isfinite(d),2),1:3);
    N = size(d,1);
    for i = 1:nq
        sel = find(qn > edges(i) & qn <= edges(i+1));
        if isempty(sel), continue, end
        if numel(sel) > ndir, sel = sel(randperm(numel(sel),ndir)); end
        ph = d*Q(sel,:).';
        Sg(i)   = Sg(i)   + sum(sum(cos(ph),1).^2 + sum(sin(ph),1).^2)/N;
        cont(i) = cont(i) + numel(sel);
    end
end
Sg = Sg./max(cont,1);
empty = cont == 0;
if any(empty)
    Sg(empty) = interp1(qg(~empty),Sg(~empty),qg(empty),'linear','extrap');
end
end
