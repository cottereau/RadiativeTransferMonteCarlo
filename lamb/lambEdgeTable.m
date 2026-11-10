function E = lambEdgeTable(vpvs)
%LAMBEDGETABLE  Rayleigh <-> 90 degree free edge coefficients vs obliquity psi.
%   E = lambEdgeTable()            homogenised medium vp/vs = 4814/2962
%   E = lambEdgeTable(vpvs)
%
%   From lambWedge25D (ppw 8, psi spread 2 deg), cached in lambEdgeTable_<vpvs>.mat next to this
%   file. Rayleigh incident on face A with angle psi between its ray and the edge normal:
%     E.psi (deg), E.T, E.R, E.P, E.S        energy fractions (T: Rayleigh on face B)
%     E.thEdges (rad), E.pP, E.pS            body energy fraction per theta bin (rows: psi),
%                                            theta from face A (0: back along A, pi/2: along B)
%   psi = 90 is added as T = 1 (grazing; 85 deg is not resolved by the 2.5-D runs).
%   Reverse path (detailed balance + reciprocity, checked at 1-5 % in 2-D): a body wave c from
%   in-plane angle theta (from face A) is converted into a Rayleigh wave on A with width
%   sigma = lambda_c pc(theta; psi)/dtheta (S: half of it per unpolarised particle), sin(psi) = a cR/c,
%   a = its cosine along the edge.
arguments
    vpvs (1,1) double = 4814/2962
end
f = fullfile(fileparts(mfilename('fullpath')), sprintf('lambEdgeTable_%.4f.mat', vpvs));
if isfile(f), E = load(f).E; return, end
ps = 0:5:80;  n = numel(ps);
E = struct('vpvs', vpvs, 'psi', [ps 90], 'T', zeros(1, n + 1), 'R', 0, 'P', 0, 'S', 0);
for i = 1:n
    o = lambWedge25D('psi', ps(i), 'ppw', 8, 'nBin', 36, 'vpvs', vpvs);
    s = o.sum;                                                      % renormalise the ~1e-4 balance
    E.T(i) = o.T/s;  E.R(i) = o.R/s;  E.P(i) = o.P/s;  E.S(i) = o.S/s;
    E.pP(i, :) = o.patP/s;  E.pS(i, :) = o.patS/s;  E.thEdges = o.thetaEdges;
end
E.T(n + 1) = 1;  E.R(n + 1) = 0;  E.P(n + 1) = 0;  E.S(n + 1) = 0;
E.pP(n + 1, :) = 0;  E.pS(n + 1, :) = 0;
b = E.P + E.S;                                                     % patterns carry P, S exactly
for i = 1:n
    E.pP(i, :) = E.pP(i, :)*E.P(i)/max(sum(E.pP(i, :)), realmin);
    E.pS(i, :) = E.pS(i, :)*E.S(i)/max(sum(E.pS(i, :)), realmin);
end
E.check = max(abs(E.T + E.R + b - 1));
save(f, 'E');
end
