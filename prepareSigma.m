function mat = prepareSigma(mat,d)

if isempty(mat.sigma)
    error(['The Differential Scattering Cross-Sections '...
        'has been not defined, please defined it usign DSCS Class'])
end

if mat.acoustics
    
    [mat.Sigma,mat.Sigmapr,mat.invcdf] = prepareSigmaOne(mat.sigma{1},d);

    % Diffusion coefficient m²/s (Eq. (5.12), Ryzhik et al, 1996)
    mat.D = mat.v^2/(d*(mat.Sigma-mat.Sigmapr));

    mat.meanFreeTime = 1/mat.Sigma;
    % the two lines below are just for homogenization of the propagation
    % code between acoustics and elastics
    mat.vp = mat.v;
    mat.vs = 0; 
else
    mat.Sigma = zeros(2);
    mat.invcdf = cell(2);

    [mat.Sigma(1,1),mat.Sigmapr(1,1),mat.invcdf{1,1}] = prepareSigmaOne(mat.sigma{1,1},d);
    [mat.Sigma(1,2),mat.Sigmapr(1,2),mat.invcdf{1,2}] = prepareSigmaOne(mat.sigma{1,2},d);
    [mat.Sigma(2,1),mat.Sigmapr(2,1),mat.invcdf{2,1}] = prepareSigmaOne(mat.sigma{2,1},d);
    [mat.Sigma(2,2),mat.Sigmapr(2,2),mat.invcdf{2,2}] = prepareSigmaOne(mat.sigma{2,2},d);
    mat.meanFreeTime = 1./sum(mat.Sigma,2);
    
    K = mat.vp/mat.vs;
    % Transport mean free paths of P & S waves
    tmfp_P = (mat.vp*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2)) + mat.vs*mat.Sigmapr(1,2) )...
                      /( (mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1))*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2))- mat.Sigmapr(1,2)*mat.Sigmapr(2,1) );
    tmfp_S = (mat.vs*(mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1)) + mat.vp*mat.Sigmapr(2,1) )...
                      /( (mat.Sigma(1,1)+mat.Sigma(1,2)-mat.Sigmapr(1,1))*(mat.Sigma(2,2)+mat.Sigma(2,1)-mat.Sigmapr(2,2))- mat.Sigmapr(1,2)*mat.Sigmapr(2,1) );
    % partial diffusion coefficients of P & S waves
    Dp = mat.vp*tmfp_P/d; 
    Ds = mat.vs*tmfp_S/d;
    % Diffusion coefficient m²/s (Eqs. (5.42) & (5.46), Ryzhik et al, 1996)
    mat.D = (Dp+2*K^3*Ds)/(1+2*K^3);

    mat.P2P = mat.Sigma(1,1)/sum(mat.Sigma(1,:),2); 
    mat.S2S = mat.Sigma(2,2)/sum(mat.Sigma(2,:),2);
end
end