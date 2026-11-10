function out = ZoeppritzFluid(theta,vp1,vs1,rho1,vp2,vs2,rho2)

% --- Material properties ---
% Material A (solid, e.g. shale)
% vp1 = 2500; % m/s
% vs1 = 1100; % m/s
% rho1 = 2200; % kg/m^3
% % Material B (fluid, e.g. water/gas)
% vp2 = 1500; % m/s (fluid velocity)
% vs2 = 0;    % <--- CHANGE: a fluid carries no S wave
% rho2 = 1000; % kg/m^3 (fluid density)
% vp2 = 343
% rho2 = 1

%theta = linspace(0, 90, 2*181);
theta_rad = theta * pi / 180;

% --- Initialise amplitude coefficients ---
Rpp = zeros(size(theta));
Rpsv = Rpp;
Tpp = Rpp;
Tpsv = Rpp;
Rsp = Rpp;
Rsvsv = Rpp;
Tsp = Rpp;
Tsvsv = Rpp;
Rsh = Rpp;
Tsh = Rpp;

E_Rpp = Rpp;
E_Rpsv = Rpp;
E_Tpp = Rpp;
E_Tpsv = Rpp;
E_P_Total = Rpp;
E_Rsp = Rpp;
E_Rsvsv = Rpp;
E_Tsp = Rpp;
E_Tsvsv = Rpp;
E_SV_Total = Rpp;
E_Rsh = Rpp;
E_Tsh = Rpp;
E_SH_Total = Rpp;

M = zeros(4,4);
N = zeros(2,2);

for iang = 1:numel(theta)

    % --- P-WAVE INCIDENCE ---
    phi0_p = theta_rad(iang);
    p = sin(phi0_p) / vp1; % ray parameter (slowness)

    phi1_p = phi0_p;
    psi1_p = asin(vs1 * p);
    phi2_p = asin(vp2 * p);
    psi2_p = asin(vs2 * p); % Isso será 0, pois vs2 = 0

    cphi1_p = cos(phi1_p);
    cpsi1_p = cos(psi1_p);
    cphi2_p = cos(phi2_p);
    cpsi2_p = cos(psi2_p); % Isso será 1

    % Matrix M (4x4) as in the solid-solid case (some terms vanish)
    M(1,1) = cphi1_p;
    M(1,2) = -sin(psi1_p);
    M(1,3) = cphi2_p;
    M(1,4) = sin(psi2_p); % 0

    M(2,1) = -sin(phi1_p);
    M(2,2) = -cpsi1_p;
    M(2,3) = sin(phi2_p);
    M(2,4) = -cpsi2_p; % -1

    M(3,1) = -cos(2*psi1_p);
    M(3,2) = vs1/vp1 * sin(2*psi1_p);
    M(3,3) = (rho2*vp2/(rho1*vp1))*cos(2*psi2_p); % cos(2*psi2_p) = 1
    M(3,4) = (rho2*vs2/(rho1*vp1))*sin(2*psi2_p); % 0

    M(4,1) = sin(2*phi1_p);
    M(4,2) = vp1/vs1 * cos(2*psi1_p); % using the corrected M(4,2)
    M(4,3) = (rho2*vp1*vs2*vs2/(rho1*vp2*vs1*vs1))*sin(2*phi2_p); % 0
    M(4,4) = -((rho2*vp1*vs2)/(rho1*vs1*vs1))*cos(2*psi2_p); % 0

    % Vector F for P incidence
    psi0_p = asin(vs1*sin(phi0_p)/vp1);
    FP(1) = cos(phi0_p);
    FP(2) = sin(phi0_p);
    FP(3) = cos(2*psi0_p);
    FP(4) = sin(2*phi0_p);

    % --- CHANGE: solve the 3x3 solid-fluid system (Rpp, Rpsv, Tpp) ---
    % Remove the 2nd row (continuity of u_x, which does not apply)
    % Remove the 4th column (Tpsv, which is 0)
    % Unknowns are [Rpp; Rpsv; Tpp]
    % Equations are (1) u_z, (3) tau_zz, (4) tau_xz = 0

    M_3x3 = [ M(1,1) M(1,2) M(1,3);
        M(3,1) M(3,2) M(3,3);
        M(4,1) M(4,2) M(4,3) ]; % M(4,3) will be 0, which is correct

    FP_3x1 = [ FP(1); FP(3); FP(4) ];

    % Resolver amplitudes P
    Ap_3x1 = M_3x3 \ FP_3x1;

    Rpp(iang) = Ap_3x1(1);
    Rpsv(iang) = Ap_3x1(2);
    Tpp(iang) = Ap_3x1(3);
    Tpsv(iang) = 0; % By definition

    % --- ENERGY COMPUTATION (P incidence) ---
    E_Rpp(iang)  = abs(Rpp(iang))^2;
    E_Rpsv(iang) = (rho1 * vs1 * real(cpsi1_p)) / (rho1 * vp1 * real(cphi1_p)) * abs(Rpsv(iang))^2;
    E_Tpp(iang)  = (rho2 * vp2 * real(cphi2_p)) / (rho1 * vp1 * real(cphi1_p)) * abs(Tpp(iang))^2;
    E_Tpsv(iang) = 0; % vs2 = 0

    if theta_rad(iang) == 0
        E_Rpsv(iang) = 0;
    end

    E_P_Total(iang) = E_Rpp(iang) + E_Rpsv(iang) + E_Tpp(iang) + E_Tpsv(iang);

    % --- SV-WAVE INCIDENCE ---
    psi0_sv = theta_rad(iang);
    p = sin(psi0_sv) / vs1;

    phi1_sv = asin(vp1 * p);
    psi1_sv = psi0_sv;
    phi2_sv = asin(vp2 * p);
    psi2_sv = asin(vs2 * p); % 0

    cphi1_sv = cos(phi1_sv);
    cpsi1_sv = cos(psi1_sv);
    cphi2_sv = cos(phi2_sv);
    cpsi2_sv = cos(psi2_sv); % 1

    % Matrix M for SV incidence
    M(1,1) = cphi1_sv;
    M(1,2) = -sin(psi1_sv);
    M(1,3) = cphi2_sv;
    M(1,4) = sin(psi2_sv); % 0

    M(2,1) = -sin(phi1_sv);
    M(2,2) = -cpsi1_sv;
    M(2,3) = sin(phi2_sv);
    M(2,4) = -cpsi2_sv; % -1

    M(3,1) = -cos(2*psi1_sv);
    M(3,2) = vs1/vp1 * sin(2*psi1_sv);
    M(3,3) = (rho2*vp2/(rho1*vp1))*cos(2*psi2_sv); % cos(2*psi2_sv) = 1
    M(3,4) = (rho2*vs2/(rho1*vp1))*sin(2*psi2_sv); % 0

    M(4,1) = sin(2*phi1_sv);
    M(4,2) = (vp1/vs1) * cos(2*psi1_sv); % <--- CORRECTION KEPT
    M(4,3) = (rho2*vp1*vs2*vs2/(rho1*vp2*vs1*vs1))*sin(2*phi2_sv); % 0
    M(4,4) = -((rho2*vp1*vs2)/(rho1*vs1*vs1))*cos(2*psi2_sv); % 0

    % Vector F for SV incidence
    FS(1) = -sin(psi0_sv);
    FS(2) = cos(psi0_sv);
    FS(3) = -vs1*sin(2*psi0_sv)/vp1;
    FS(4) = vp1*cos(2*psi0_sv)/vs1;

    % --- CHANGE: solve the 3x3 solid-fluid system (Rsp, Rsvsv, Tsp) ---
    % Remove the 2nd row (continuity of u_x) and the 4th column (Tsvsv)
    M_3x3 = [ M(1,1) M(1,2) M(1,3);
        M(3,1) M(3,2) M(3,3);
        M(4,1) M(4,2) M(4,3) ]; % M(4,3) será 0

    FS_3x1 = [ FS(1); FS(3); FS(4) ];

    % Resolver amplitudes SV
    As_3x1 = M_3x3 \ FS_3x1;

    Rsp(iang) = As_3x1(1);
    Rsvsv(iang) = As_3x1(2);
    Tsp(iang) = As_3x1(3);
    Tsvsv(iang) = 0; % By definition

    % --- ENERGY COMPUTATION (SV incidence) ---
    E_Rsp(iang)   = (rho1 * vp1 * real(cphi1_sv)) / (rho1 * vs1 * real(cpsi1_sv)) * abs(Rsp(iang))^2;
    E_Rsvsv(iang) = abs(Rsvsv(iang))^2;
    E_Tsp(iang)   = (rho2 * vp2 * real(cphi2_sv)) / (rho1 * vs1 * real(cpsi1_sv)) * abs(Tsp(iang))^2;
    E_Tsvsv(iang) = 0; % vs2 = 0

    if theta_rad(iang) == 0, E_Rsp(iang) = 0; end

    E_SV_Total(iang) = E_Rsp(iang) + E_Rsvsv(iang) + E_Tsp(iang) + E_Tsvsv(iang);

    % --- SH-WAVE INCIDENCE ---
    psi0_sh = theta_rad(iang);
    p = sin(psi0_sh) / vs1;

    psi1_sh = psi0_sh;
    psi2_sh = asin(vs2 * p); % 0

    cpsi1_sh = cos(psi1_sh);
    cpsi2_sh = cos(psi2_sh); % 1

    % --- CHANGE: solid-fluid for SH ---
    % The fluid carries no SH wave (vs2 = 0).
    % The boundary condition is that the shear traction (tau_yz)
    % vanishes at the interface.
    % tau_yz(inc) + tau_yz(ref) = 0
    % Z1 * (1 - Rsh) = 0, where Z1 = rho1*vs1*cpsi1_sh
    % This forces Rsh = 1 (total reflection) and Tsh = 0 (no transmission).

    Rsh(iang) = 1.0;
    Tsh(iang) = 0.0;

    % Handle the division by zero at 90 deg (cpsi1_sh = 0) in the energy
    if abs(cpsi1_sh) < 1e-10
        E_Rsh(iang) = 1.0;
        E_Tsh(iang) = 0.0;
    else
        % --- ENERGY COMPUTATION (SH incidence) ---
        E_Rsh(iang) = abs(Rsh(iang))^2;
        E_Tsh(iang) = (rho2 * vs2 * real(cpsi2_sh)) / (rho1 * vs1 * real(cpsi1_sh)) * abs(Tsh(iang))^2; % Será 0
    end

    E_SH_Total(iang) = E_Rsh(iang) + E_Tsh(iang);
end

out.E_Rpp   = E_Rpp;
out.E_Rpsv  = E_Rpsv;
out.E_Tpp   = E_Tpp;
out.E_Tpsv  = E_Tpsv;
out.E_Rsvsv = E_Rsvsv; 
out.E_Rsp   = E_Rsp;
out.E_Tsvsv = E_Tsvsv; 
out.E_Tsp   = E_Tsp;
out.E_Rsh   = E_Rsh;
out.E_Tsh   = E_Tsh;

out.A_Rpp   = Rpp;
out.A_Rpsv  = Rpsv;
out.A_Tpp   = Tpp;
out.A_Tpsv  = Tpsv;
out.A_Rsvsv = Rsvsv; 
out.A_Rsp   = Rsp;
out.A_Tsvsv = Tsvsv; 
out.A_Tsp   = Tsp;
out.A_Rsh   = Rsh;
out.A_Tsh   = Tsh;

if 0
    % --- PLOTS ---
    % (Amplitude and energy plots kept as in the original)
    % Zeroed coefficients (Tpsv, Tsvsv, Tsh) appear as a line at 0.
    % Coefficient Rsh appears as a line at 1.

    figure('Name', 'SH');
    plot(theta, abs(Rsh), 'b', 'LineWidth', 1.5);
    hold on;
    plot(theta, abs(Tsh), 'r', 'LineWidth', 1.5);
    title('Amplitudes (SH-Incidence) - Solid-Fluid');
    xlabel('Incidence angle (deg)');
    ylabel('Magnitude da Amplitude');
    legend('Rsh', 'Tsh');
    grid on;
    hold off;

    figure('Name', 'SV');
    plot(theta, abs(Rsvsv), 'b', 'LineWidth', 1.5);
    hold on;
    plot(theta, abs(Rsp), 'r', 'LineWidth', 1.5);
    plot(theta, abs(Tsvsv), 'g', 'LineWidth', 1.5);
    plot(theta, abs(Tsp), 'm', 'LineWidth', 1.5);
    title('Amplitudes (SV-Incidence) - Solid-Fluid');
    xlabel('Incidence angle (deg)');
    ylabel('Magnitude da Amplitude');
    legend('Rsvsv', 'Rsp', 'Tsvsv', 'Tsp');
    grid on;
    hold off;

    figure('Name', 'P');
    plot(theta, abs(Rpp), 'b', 'LineWidth', 1.5);
    hold on;
    plot(theta, abs(Rpsv), 'r', 'LineWidth', 1.5);
    plot(theta, abs(Tpp), 'g', 'LineWidth', 1.5);
    plot(theta, abs(Tpsv), 'm', 'LineWidth', 1.5);
    title('Amplitudes (P-Incidence) - Solid-Fluid');
    xlabel('Incidence angle (deg)');
    ylabel('Magnitude da Amplitude');
    legend('Rpp', 'Rpsv', 'Tpp', 'Tpsv');
    grid on;
    hold off;

    
    % --- NOVOS PLOTS DE ENERGIA ---
    figure('Name', 'Energy (P incidence)');
    plot(theta, E_Rpp, 'b', 'LineWidth', 1.5);
    hold on;
    plot(theta, E_Rpsv, 'r', 'LineWidth', 1.5);
    plot(theta, E_Tpp, 'g', 'LineWidth', 1.5);
    plot(theta, E_Tpsv, 'm', 'LineWidth', 1.5);
    plot(theta, E_P_Total, 'k--', 'LineWidth', 2.5); % Soma total
    title('Energy partitioning (P incidence) - Solid-Fluid');
    xlabel('Incidence angle (deg)');
    ylabel('Energy coefficient');
    legend('E_{Rpp}', 'E_{Rpsv}', 'E_{Tpp}', 'E_{Tpsv}', 'E_{Total}');
    ylim([-0.1 1.1]);
    grid on;
    hold off;

    figure('Name', 'Energy (SV incidence)');
    plot(theta, E_Rsvsv, 'b', 'LineWidth', 1.5);
    hold on;
    plot(theta, E_Rsp, 'r', 'LineWidth', 1.5);
    plot(theta, E_Tsvsv, 'g', 'LineWidth', 1.5);
    plot(theta, E_Tsp, 'm', 'LineWidth', 1.5);
    plot(theta, E_SV_Total, 'k--', 'LineWidth', 2.5); % Soma total
    title('Energy partitioning (SV incidence) - Solid-Fluid');
    xlabel('Incidence angle (deg)');
    ylabel('Energy coefficient');
    legend('E_{Rsvsv}', 'E_{Rsp}', 'E_{Tsvsv}', 'E_{Tsp}', 'E_{Total}');
    ylim([-0.1 1.1]);
    grid on;
    hold off;

    figure('Name', 'Energy (SH incidence)');
    plot(theta, E_Rsh, 'b', 'LineWidth', 1.5);
    hold on;
    plot(theta, E_Tsh, 'r', 'LineWidth', 1.5);
    plot(theta, E_SH_Total, 'k--', 'LineWidth', 2.5); % Soma total
    title('Energy partitioning (SH incidence) - Solid-Fluid');
    xlabel('Incidence angle (deg)');
    ylabel('Energy coefficient');
    legend('E_{Rsh}', 'E_{Tsh}', 'E_{Total}');
    ylim([-0.1 1.1]);
    grid on;
    hold off;
end