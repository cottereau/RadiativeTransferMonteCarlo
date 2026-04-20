function [mean_lam, std_lam, mean_mu, std_mu, mean_rho, std_rho] = ...
    convert_kmr_to_lmr(mean_kappa, std_kappa, mean_mu_in, std_mu_in, mean_rho_in, std_rho_in, corr_K_mu)
%
%   converts the mean values and standard deviations of the isotropic
%   material parameters (\kappa, \mu, \rho) into the corresponding
%   statistics for (\lambda, \mu, \rho), using 3D isotropic relation:
%
%       \lambda = \kappa - (2/3)\mu
%
%   INPUTS
%       mean_kappa   - Mean bulk modulus \kappa
%       std_kappa    - Standard deviation of \kappa
%       mean_mu_in   - Mean shear modulus \mu
%       std_mu_in    - Standard deviation of \mu
%       mean_rho_in  - Mean density \rho
%       std_rho_in   - Standard deviation of \rho
%       corr_K_mu    - Correlation coefficient between \kappa and \mu
%
%   OUTPUTS
%       mean_lam     - Mean Lamé parameter \lambda
%       std_lam      - Standard deviation of \lambda
%       mean_mu      - Mean shear modulus \mu (unchanged)
%       std_mu       - Standard deviation of \mu (unchanged)
%       mean_rho     - Mean density \rho (unchanged)
%       std_rho      - Standard deviation of \rho (unchanged)
%
%   NOTES
%       The transformation is based on linear isotropic elasticity in 3D.
%       Since \lambda is a linear combination of \kappa and \mu, the mean
%       value is obtained exactly as
%
%           E[\lambda] = E[\kappa] - (2/3)E[\mu].
%
%       The variance is computed as
%
%           Var[\lambda] = Var[\kappa] + (4/9)Var[\mu]
%                        - (4/3)Cov(\kappa,\mu).
%
%   WARNING
%       If \kappa and \mu are correlated, the returned value of std_lam
%       will generally be incorrect unless the covariance term is included.
%
%   EXAMPLE
%       [mean_lam, std_lam, mean_mu, std_mu, mean_rho, std_rho] = ...
%           convert_kmr_to_lmr(50e9, 2e9, 30e9, 1e9, 7800, 50);
%

% default value for the covariance of kappa and mu is 0
if nargin < 7
    corr_K_mu = 0;
    warning('Assuming kappa and mu are uncorrelated: corr(kappa,mu)=0');
end

mean_mu  = mean_mu_in;
std_mu   = std_mu_in;
mean_rho = mean_rho_in;
std_rho  = std_rho_in;

mean_lam = mean_kappa - (2/3)*mean_mu;
var_lam  = std_kappa^2 + (4/9)*std_mu^2 - (4/3)*std_kappa*std_mu*corr_K_mu;
std_lam = sqrt(var_lam);

if var_lam < 0
    error('Computed variance of lambda is negative. Check inputs.');
end

if mean_kappa <= 0 || mean_mu <= 0 || mean_rho <= 0
    warning('Nonpositive mean values detected in kappa, mu, or rho.');
end

if 3*mean_lam + 2*mean_mu <= 0
    warning('Elastic stability condition 3*lambda + 2*mu > 0 is violated.');
end

end