function mat = prepareSigma(mat,d)
if isfield(mat,'sigma')
    [mat.Sigma,mat.invcdf] = prepareSigmaOne(mat.sigma,d);
    mat.meanFreeTime = 1/mat.v/mat.Sigma;
end
if isfield(mat,'sigmaPP')
    [mat.SigmaPP,mat.invcdfPP] = prepareSigmaOne(mat.sigmaPP,d);
end
if isfield(mat,'sigmaPS')
    [mat.SigmaPS,mat.invcdfPS] = prepareSigmaOne(mat.sigmaPS,d);
    mat.meanFreeTimeP = 1/mat.vp/(mat.SigmaPP + mat.SigmaPS); % This should be changed (mean free time => to be multiplied by vp)
end
if isfield(mat,'sigmaSS')
    [mat.SigmaSS,mat.invcdfSS] = prepareSigmaOne(mat.sigmaSS,d);
end
if isfield(mat,'sigmaSP')
    [mat.SigmaSP,mat.invcdfSP] = prepareSigmaOne(mat.sigmaSP,d);
    mat.meanFreeTimeS = 1/mat.vs/(mat.SigmaSS + mat.SigmaSP); % This should be changed (mean free time => to be multiplied by vs)
    if 2*(mat.vp/mat.vs)^2*mat.SigmaSP~=mat.SigmaPS    % This should be changed in 3D (2K^3)
        disp('error in the Sigma SP/PS relation:')
        disp(['SigmaPS = ' num2str(mat.SigmaSP) ])
        disp(['2(vP/vS)^2 SigmaSP = ' num2str(2*(mat.vp/mat.vs)^2*mat.SigmaSP) ])
    end
    % Question: why do you need to calculate these probabilities? We rather
    % need mat.P2S and mat.S2P no?
    mat.P2P = mat.SigmaPP/(mat.SigmaPP+mat.SigmaPS); % mat.SigmaPS should it be rather mat.SigmaSP? In Ryzhik/Baydoun Sigma_ij is j to i as oposed to in Margerin/Turner/Weaver
    mat.S2S = mat.SigmaSS/(mat.SigmaSS+mat.SigmaSP); % mat.SigmaSP should it be rather mat.SigmaPS? The same comment ...
end
end


