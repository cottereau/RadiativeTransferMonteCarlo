function [k, w, cells] = windowNodes(breaks, n, peaks)
%WINDOWNODES  Quadrature nodes on [0, k_m] for integrands with branch-point and
%   near-pole singularities.
%   [k, w] = lamb.quad.windowNodes(breaks, n)
%   [k, w] = lamb.quad.windowNodes(breaks, n, peaks)
%
%   breaks : branch points [0, k_1, ..., k_m]; the window is split at every one of them.
%            On a segment [ka, kb] touching a branch point the substitution
%            k^2 = ka^2 + (kb^2 - ka^2) sin^2 t (k = kb sin t if ka = 0) removes the
%            1/sqrt singularities at both ends (eq:subst); Gauss-Legendre in t.
%   n      : Gauss nodes per sub-segment.
%   peaks  : P-by-2 [a b] Lorentzian peaks (near-real leaky poles a + 1i b). A segment that
%            contains a is further split at a -+ b*10.^(0:J) and those interior pieces use
%            plain Gauss-Legendre in k, which resolves peaks of any width b. A peak with a
%            at a segment end (a branch point) grades that end (lamb.quad.peaksFor).
%   k, w   : column vectors; int f(k) dk ~ sum(w .* f(k)).
%   cells  : N-by-5 [tLo tHi ka kb type], the quadrature cell of each node in its segment's
%            mapped variable t: edges at the cumulative Gauss weights, so cell n has width
%            w_t(n) in t and contains node n. type 1: k^2 = ka^2 + (kb^2-ka^2) sin^2 t;
%            type 2: k = (ka+kb)/2 + (kb-ka)/2 t. Used to sample k inside the node's own
%            share of the quadrature mass (lamb.quad.cellSample).
if nargin < 3 || isempty(peaks), peaks = zeros(0, 2); end
breaks = unique(breaks(:)).';
[t, wt] = gaussLegendre(n, 0, pi/2);
[x, wx] = gaussLegendre(n, -1, 1);
k = [];  w = [];  cells = zeros(0, 5);
eT = [0; cumsum(wt)];  eX = [-1; -1 + cumsum(wx)];            % cumulative-weight edges
for s = 1:numel(breaks) - 1
    lo = breaks(s);  hi = breaks(s+1);
    pts = [lo, hi];
    for p = find(peaks(:,1) > lo & peaks(:,1) <= hi).'          % a == hi: grade toward a branch point
        a = peaks(p,1);  b = peaks(p,2);
        g = b*10.^(0:ceil(log10((hi - lo)/b)));
        pts = [pts, a - g, a + g];                                  %#ok<AGROW>
    end
    pts = unique(pts(pts >= lo & pts <= hi));
    for j = 1:numel(pts) - 1
        ka = pts(j);  kb = pts(j+1);
        if j == 1 || j == numel(pts) - 1                            % touches a branch point
            if ka == 0
                kk = kb*sin(t);  dk = kb*cos(t);
            else
                kk = sqrt(ka^2 + (kb^2 - ka^2)*sin(t).^2);
                dk = (kb^2 - ka^2)*sin(t).*cos(t)./kk;
            end
            ww = dk.*wt;
            cc = [eT(1:end-1), eT(2:end), repmat([ka kb 1], n, 1)];
        else                                                        % interior (peak) piece
            kk = (ka + kb)/2 + (kb - ka)/2*x;
            ww = (kb - ka)/2*wx;
            cc = [eX(1:end-1), eX(2:end), repmat([ka kb 2], n, 1)];
        end
        k = [k; kk];  w = [w; ww];  cells = [cells; cc];            %#ok<AGROW>
    end
end
end

function [x, w] = gaussLegendre(n, a, b)
% Golub-Welsch nodes and weights on [a, b].
i = (1:n-1)';
beta = i./sqrt(4*i.^2 - 1);
[V, D] = eig(diag(beta, 1) + diag(beta, -1));
[x, idx] = sort(diag(D));
w = 2*V(1, idx)'.^2;
x = (b - a)/2*x + (a + b)/2;
w = (b - a)/2*w;
end
