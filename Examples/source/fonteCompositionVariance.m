function budget = fonteCompositionVariance(bandEnergies,unitEstimatorCovariance)
%FONTECOMPOSITIONVARIANCE Variance of a physical spectral-energy sum.
% bandEnergies(j,beta) are weighted physical energies in J. Covariance is
% for the matching unit-response estimators, ordered by MATLAB (:), and
% already includes their particle/batch counts. No independence is assumed.
validateattributes(bandEnergies,{'numeric'},{'real','finite','2d','nonnegative'});
if size(bandEnergies,2)~=2 || isempty(bandEnergies)
    error('Fonte:VarianceInput','Use an N-by-2 [P S] band-energy array.');
end
a=double(bandEnergies(:)); n=numel(a);
validateattributes(unitEstimatorCovariance,{'numeric'}, ...
    {'real','finite','size',[n n]});
C=double(unitEstimatorCovariance);
if norm(C-C','fro')>1e-10*max(1,norm(C,'fro'))
    error('Fonte:VarianceInput','Estimator covariance must be symmetric.');
end
C=(C+C')/2;
if min(eig(C)) < -1e-10*max(1,norm(C,2))
    error('Fonte:VarianceInput','Estimator covariance must be positive semidefinite.');
end
diagonal=a.^2.*diag(C);
variance=max(0,a'*C*a);
budget=struct('variance',variance,'standardError',sqrt(variance), ...
    'diagonalContributions',reshape(diagonal,size(bandEnergies)), ...
    'covarianceContribution',variance-sum(diagonal), ...
    'energyUnits','J','varianceUnits','J^2');
end
