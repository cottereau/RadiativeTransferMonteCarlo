function [med_lam, std_lam, med_mu, std_mu, med_rho, std_rho] = ...
    convert_kmr_para_lmr(med_k, std_k, med_mu_in, std_mu_in, med_rho_in, std_rho_in)
% CONVERTE_KMR_PARA_LMR
% Converts K (Bulk), Mu, Rho -> Lambda, Mu, Rho
% with uncertainty propagation assuming statistical independence.

warning('Assuming kappa and mu are uncorrelated: cov(kappa,mu)=0')
cov_K_mu = 0;
% Mu and Rho are copied directly
med_mu  = med_mu_in;
std_mu  = std_mu_in;
med_rho = med_rho_in;
std_rho = std_rho_in;

% ====================== Mean Calculation ======================
med_lam = med_k - (2/3) * med_mu;

% ====================== Uncertainty Propagation ======================
% lambda = K - (2/3)mu -> partial derivative with respect to mu = -2/3
var_lam = std_k^2 + (4/9) * std_mu^2 - 2*(2/3)*cov_K_mu;
std_lam = sqrt(var_lam);

% ====================== Validations ======================
if med_k < 0 || med_mu < 0 || med_rho <= 0
    warning('Negative values detected in K or Mu, or invalid density.');
end

if med_lam < -1e-10*abs(med_k)  % allows slightly negative lambda
    warning('Resulting lambda is negative.');
end

end