function mat = prepareSigma(mat,d)
if mat.acoustics
    [mat.Sigma,mat.invcdf] = prepareSigmaOne(mat.sigma,d);
    mat.meanFreeTime = 1/mat.v/mat.Sigma;
    % the two lines below are just for homogenization of the propagation
    % code between acoustics and elastics
    mat.vp = mat.v;
    mat.vs = 0; 
else
    mat.Sigma = zeros(2);
    mat.invcdf = cell(2);
    [mat.Sigma(1,1),mat.invcdf{1,1}] = prepareSigmaOne(mat.sigma{1,1},d);
    [mat.Sigma(1,2),mat.invcdf{1,2}] = prepareSigmaOne(mat.sigma{1,2},d);
    [mat.Sigma(2,1),mat.invcdf{2,1}] = prepareSigmaOne(mat.sigma{2,1},d);
    [mat.Sigma(2,2),mat.invcdf{2,2}] = prepareSigmaOne(mat.sigma{2,2},d);
    mat.meanFreeTime = 1./[mat.vp;math.vs]./sum(mat.Sigma,2);
    % Question: why do you need to calculate these probabilities? We rather
    % need mat.P2S and mat.S2P no?
    mat.P2P = mat.Sigma(1,1)/sum(mat.Sigma,2); % mat.SigmaPS should it be rather mat.SigmaSP? In Ryzhik/Baydoun Sigma_ij is j to i as oposed to in Margerin/Turner/Weaver
    mat.S2S = mat.Sigma(2,2)/sum(mat.Sigma,2); % mat.SigmaSP should it be rather mat.SigmaPS? The same comment ...
end
end


