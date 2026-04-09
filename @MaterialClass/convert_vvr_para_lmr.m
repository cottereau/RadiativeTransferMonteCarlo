function [med_lam, std_lam, med_mu, std_mu, med_rho, std_rho] = ...
    convert_vvr_para_lmr(med_vp, std_vp, med_vs, std_vs, med_rho_in, std_rho_in)
% CONVERTE_VVR_PARA_LMR
% Converts Vp, Vs, Rho -> Lambda, Mu, Rho with uncertainty propagation
% Uses the Delta Method

% Rho remains unchanged
med_rho = med_rho_in;
std_rho = std_rho_in;

% ====================== Mean Calculation ======================
med_mu  = med_rho * med_vs^2;
med_lam = med_rho * (med_vp^2 - 2 * med_vs^2);

% ====================== Uncertainty Propagation ======================

% --- For Mu = rho * vs^2 ---
dmu_drho = med_vs^2;
dmu_dvs  = 2 * med_rho * med_vs;

var_mu   = (dmu_drho^2 * std_rho^2) + (dmu_dvs^2 * std_vs^2);
std_mu   = sqrt(var_mu);

% --- For Lambda = rho*(vp^2 - 2*vs^2) ---
dlam_drho = (med_vp^2 - 2*med_vs^2);
dlam_dvp  = 2 * med_rho * med_vp;
dlam_dvs  = -4 * med_rho * med_vs;   % derivative with respect to vs

var_lam   = (dlam_drho^2 * std_rho^2) + ...
    (dlam_dvp^2  * std_vp^2)  + ...
    (dlam_dvs^2  * std_vs^2);
std_lam   = sqrt(var_lam);

% ====================== Validations and Warnings ======================
if med_vs <= 0 || med_vp <= 0 || med_rho <= 0
    warning('Invalid physical values: velocities or density cannot be negative or zero.');
end
if med_vp < med_vs
    warning('Vp less than Vs detected. This is physically impossible in isotropic materials.');
end

end