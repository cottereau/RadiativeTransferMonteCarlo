function P = scatterParticle(acoustics,mat,P)
%function [dtot,p,v,meanFreePath] = scatterParticle(d,mat,p,v,meanFreePath)
N = P.N;
rd = rand(N,1);
% acoustics
if acoustics
    th = mat.invcdf(rd);
% elastics
else
    p0 = p;
    % change polarization
    probabilityOfChange = mat.P2P*p + mat.S2S*(1-p);
    change = rand(N,1)>probabilityOfChange;
    p(change)=~p(change);
    % velocity for each particle
    v(change&p) = mat.vp;
    v(change&~p) = mat.vs;
    % meanFreePath for each particle
    meanFreePath(change&p) = mat.meanFreePathP;
    meanFreePath(change&~p) = mat.meanFreePathS;    
    % change direction depending on polarizations
    th = zeros(N,1);
    th(p0&p) = mat.invcdfPP(rand(nnz(p0&p),1));
    th(p0&~p) = mat.invcdfPS(rand(nnz(p0&~p),1));
    th(~p0&p) = mat.invcdfSP(rand(nnz(~p0&p),1));
    th(~p0&~p) = mat.invcdfSS(rand(nnz(~p0&~p),1));
end
%P.d = P.d+th;
P.costheta = cos(acos(P.costheta)+th);

