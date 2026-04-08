close all
clearvars
clc

titlecase = '3D acoustic case with isotropic scattering';
disp(['Testing ' titlecase ' ...']);

% input data
geometry = struct( 'dimension', 3 );

source = struct( 'numberParticles', 1e6, ...
                 'position', [0 0 0], ...
                 'lambda', 2e-4 );

material = MaterialClass.preset(1);

observation = struct( 'x', 0:0.1:10, ...
                      'y', [-pi pi], ...
                      'z', [-pi/2 pi/2], ...
                      'directions', [0 pi], ...
                      'time', 0:0.05:20);

inds = [20 50 80]; % index of the desired observation points

% running our code, Monte Carlo-based
obs = radiativeTransfer( geometry, source, material, observation );

Eus = squeeze(obs.energyDensity);

% running Yoshimoto's Monte Carlo-based approach
EY = Comparison.randomWalkYoshimoto( geometry, source, material, observation );

% computing Paasschens solution
[EP, Ediff] = Comparison.analyticalPaasschens( material, observation, geometry );

% visual comparison
figure; hold on; grid on; box on;
h1 = semilogy( obs.t, Eus(inds,:), '-k' );
h2 = semilogy( obs.t, EY(inds,:), '-r' );
h3 = semilogy( obs.t, EP(inds,:), '-b' );
h4 = semilogy( obs.t, Ediff(inds,:), ':r' );
set(gca, 'YScale', 'log'); ylim([1e-5 1]);
legend( [h1(1), h2(1), h3(1), h4(1)], ...
    {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
    'Analytical (Paasschens, 1997)', 'Diffusion approximation'}, ...
    'FontSize',12);
xlabel('Lapse Time [s]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);
%%
inds = [60 150 240]; 
% visual comparison
figure; hold on; grid on; box on;
h1 = semilogy( obs.x, Eus(:,inds), '-k' );
h2 = semilogy( obs.x, EY(:,inds), '-r' );
h3 = semilogy( obs.x, EP(:,inds), '-b' );
h4 = semilogy( obs.x, Ediff(:,inds), ':r' );
set(gca, 'YScale', 'log'); ylim([1e-5 1]);
legend( [h1(1), h2(1), h3(1), h4(1)], ...
    {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
    'Analytical (Paasschens, 1997)', 'Diffusion approximation'}, ...
    'FontSize',12);
xlabel('Distance [m]');
ylabel('Integrated energy density at different source-station distances')
title(titlecase);
%%

mu_a = 0;
speed = 1;
D = material.Diffusivity;
rr = obs.x;
t = obs.t;
pdefun_handle = @(x, t_val, u, dudx) pdefun(x, t_val, u, dudx, D, speed, mu_a);

CI = squeeze(obs.energyDensity(:,1,2));
icfun_handle = @(x) icfun(x, rr, CI);
sol = pdepe(2, pdefun_handle, icfun_handle, @bcfun, rr, t);
sol = sol'/(pi/2);
h5 = semilogy( obs.x, sol(:,inds+1), ':k' );
legend( [h1(1), h2(1), h3(1), h4(1), h5(1)], ...
    {'Monte Carlo (our code)', 'Monte Carlo (Yoshimoto 2000)', ...
    'Analytical (Paasschens, 1997)', 'Diffusion approximation','Diffusion calc'}, ...
    'FontSize',12);

function [c, f, s] = pdefun(~, ~, u, dudx, D, speed, mu_a)
c = 1;  % Coeficiente temporal: 1, como explicado
f = speed * D * dudx;  % Fluxo escalado pela velocidade (correção principal)
s = -speed * mu_a * u;  % Se mu_a > 0, inclui absorção; senão s=0
end

function u0 = icfun(x, r_vec, RTT_t0)
u0 = interp1(r_vec, RTT_t0, x, 'linear', 0);
end

function [pl, ql, pr, qr] = bcfun(~, ul, ~, ur, ~)
% xL = borda esquerda (r=0), xR = borda direita (r=max)
% Em r=0 (Simetria): Fluxo deve ser zero (dE/dr = 0)
pl = ul - 0;  % Pequena clareza: para Dirichlet se precisar, mas aqui é Neumann f=0
ql = 1;       % 0 + 1*f = 0 => f=0 (fluxo nulo)

% Em r=R_max: Assumindo fronteira aberta (energia zero longe, Dirichlet u=0)
pr = ur;      % ur + 0*f = 0 => ur=0
qr = 0;
% Alternativa: se quiser Neumann (fluxo zero para domínio infinito), use pr=0, qr=1 => 0 + 1*f=0 => f=0
end