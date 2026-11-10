%% Water layer over a monodisperse sphere medium, P source at the water surface
%
%  Geometry (z positive downwards):
%
%        z = 0     ---- free surface (water / air) --------------  * source
%                  water, no volume scattering
%        z = H     ==== water / composite interface ============
%                  elastic composite: monodisperse spheres,
%                  Percus-Yevick correlation, volume scattering
%        z -> inf  (half-space, no bottom boundary)
%
%  The water is a genuine fluid (vs = 0), so it carries P only, and the
%  interface is handled by MaterialInterfaceSolidFluid: a P wave arriving
%  from the water converts into P and SV in the composite, and shear energy
%  scattered back up cannot re-enter the water as shear.
%
%  The free surface at z = 0 is a fluid-vacuum boundary.
%  propagateParticleMultiMaterial reflects it specularly and never calls
%  Zoeppritz there -- a fluid against vacuum makes the Zoeppritz boundary
%  matrix singular (two zero-shear media at once), and there is nothing else
%  a fluid free surface can do anyway: no S to convert into.
%
%  The composite scattering comes from the MONODISPERSE correlation:
%  SpectralLaw = 'monodispersesphere' solves the Percus-Yevick
%  Ornstein-Zernike closure for hard spheres and builds the PSDF from the
%  resulting S(k) (Torquato & Stell 1985, doi:10.1063/1.448475). That is a real packing
%  correlation, not independent scatterers: S(k) suppresses forward
%  scattering at the packing fraction used here.
%
%  See also PROPAGATEPARTICLEMULTIMATERIAL, MATERIALINTERFACESOLIDFLUID,
%  MATERIALCLASS/MONODISPERSESPHERE.

clear; clc; rng(7);

%% 1. Parameters
f0    = 500e3;          % Hz
H     = 0.040;          % m, water layer thickness
Dsph  = 8e-3;           % m, sphere diameter
phi   = 0.35;           % volume fraction of spheres
zsrc  = 0.002;          % m, source depth below the free surface
N     = 1e6/4;          % particles
Tend  = 200e-6;         % s
nOut  = 150;            % output times

%% 2. Water (fluid, no volume scattering)
water = MaterialClass;
water.d = 3;  water.acoustics = false;
water.vp = 1480;  water.vs = 0;  water.rho = 1000;
water.Frequency = f0;
water.timeSteps = 0;
water.meanFreeTime = [Inf; Inf];       % the layer itself does not scatter
water.P2P = 1;  water.S2S = 1;
zeroInv = @(u) zeros(size(u));
water.invcdf = {zeroInv zeroInv; zeroInv zeroInv};

%% 3. Composite: monodisperse spheres, Percus-Yevick correlation
% Effective velocities for transport; the fluctuation statistics of the
% two-phase medium drive the scattering.
comp = MaterialClass;
comp.d = 3;  
comp.acoustics = false;
comp.vp = 1000; 
comp.vs = comp.vp/sqrt(3);  
comp.rho = 700;
comp.Frequency = f0;
comp.coefficients_of_variation = [0.42 0.50 0.089];   % CV of lambda, mu, rho
comp.correlation_coefficients  = [1 1 1];             % one indicator drives all three
comp.SpectralLaw  = 'monodispersesphere';
comp.SpectralParam = struct('rhoS', phi, 'Diam', Dsph);
comp.getPSDF;                                          % Percus-Yevick -> PSDF
comp = MaterialClass.prepareSigma(comp, 3);

matArray = {water, comp};

fprintf('\n=== Water over monodisperse sphere medium ===\n');
fprintf('f0 = %.0f kHz\n', f0/1e3);
fprintf('water : vp = %.0f m/s, layer H = %.1f mm, lambda_P = %.2f mm\n', ...
    water.vp, 1e3*H, 1e3*water.vp/f0);
fprintf('spheres: D = %.1f mm, phi = %.2f, k_P*D = %.2f (resonant regime)\n', ...
    1e3*Dsph, phi, 2*pi*f0/comp.vp*Dsph);
fprintf('composite: vp = %.0f, vs = %.0f m/s, Lc = %.2f mm\n', ...
    comp.vp, comp.vs, 1e3*comp.CorrelationLength);
fprintf('  mean free path  P = %.1f mm, S = %.1f mm\n', ...
    1e3*comp.vp*comp.meanFreeTime(1), 1e3*comp.vs*comp.meanFreeTime(2));
fprintf('  P2P = %.3f, S2S = %.3f\n', comp.P2P, comp.S2S);

%% 4. Geometry
geometry = struct('dimension',3,'frame','cartesian');
geometry.bnd(1) = struct('dir',3,'val',0,'type','reflective','mat',[]);  % water / air
geometry.bnd(2) = struct('dir',3,'val',H,'type','','mat',[1 2]);         % water / composite

%% 5. Source: isotropic P point source near the top of the water
% The upgoing half reflects off the free surface almost immediately, which
% is what a source just under a water surface really does.
P.N = N;
P.x = [zeros(N,2), zsrc*ones(N,1)];
c = 2*rand(N,1)-1;  ph = 2*pi*rand(N,1);  s = sqrt(1-c.^2);
P.dir  = [s.*cos(ph), s.*sin(ph), c];
P.perp = [-sin(ph), cos(ph), zeros(N,1)];
P.p    = true(N,1);                    % P source
P.alive = true(N,1);
P.t = zeros(N,1);
P.matIdx = ones(N,1);                  % starts in the water

%% 6. Propagate, recording energy by region and mode
tVec = linspace(0, Tend, nOut);
Ebud = zeros(nOut,3);                  % water P | composite P | composite S

zEdges = linspace(0, 0.30, 121);  zc = 0.5*(zEdges(1:end-1)+zEdges(2:end));
Ez = zeros(nOut, numel(zc), 2);        % depth profile, BOTH media: P, S

% Receiver array on the surface: 128 sensors, 64 each side of the source, in
% the y = 0 plane. They sit inside the water, so they can only ever record P
% -- which is the point: everything they see after the direct arrival came
% back up out of the heterogeneous medium.
dxs    = 1e-3;                                   % sensor pitch, m
xs     = dxs*[-(63.5:-1:0.5), (0.5:1:63.5)];     % 128 positions, none at x = 0
sEdges = [xs - dxs/2, xs(end)+dxs/2];
ysl    = 5e-3;                                   % half-thickness in y of a sensor cell
zsens  = 5e-3;                                   % sensor cells span 0 <= z <= zsens
Esens  = zeros(nOut, numel(xs));

xEdges = linspace(-0.15, 0.15, 121);
zgEdges = linspace(0, 0.30, 121);
xg = 0.5*(xEdges(1:end-1)+xEdges(2:end));
zg = 0.5*(zgEdges(1:end-1)+zgEdges(2:end));
slab = 0.01;
frames = zeros(nOut, numel(xg), numel(zg), 2);

fprintf('\npropagating %g particles over %d output times...\n', N, nOut);
tic
for it = 1:nOut
    if it > 1
        P = propagateParticleMultiMaterial(matArray, geometry, P, tVec(it));
    end
    inW = P.matIdx==1 & P.alive;
    inC = P.matIdx==2 & P.alive;
    Ebud(it,:) = [nnz(inW), nnz(inC & P.p), nnz(inC & ~P.p)]/N;

    ib = discretize(P.x(:,3), zEdges);
    okP = P.alive &  P.p & ~isnan(ib);
    okS = P.alive & ~P.p & ~isnan(ib);
    Ez(it,:,1) = accumarray(ib(okP), 1, [numel(zc) 1])/N;
    Ez(it,:,2) = accumarray(ib(okS), 1, [numel(zc) 1])/N;

    % Receiver array
    inCell = P.alive & abs(P.x(:,2)) < ysl & P.x(:,3) >= 0 & P.x(:,3) <= zsens;
    is = discretize(P.x(:,1), sEdges);
    hit = inCell & ~isnan(is);
    if any(hit)
        Esens(it,:) = accumarray(is(hit), 1, [numel(xs) 1]).'/N;
    end

    sl = abs(P.x(:,2)) < slab & P.alive;
    frames(it,:,:,1) = histcounts2(P.x(sl& P.p,1), P.x(sl& P.p,3), xEdges, zgEdges)/N;
    frames(it,:,:,2) = histcounts2(P.x(sl&~P.p,1), P.x(sl&~P.p,3), xEdges, zgEdges)/N;

    if mod(it, max(1,round(nOut/10)))==0
        fprintf('  %3.0f%%  water %.3f | comp P %.3f | comp S %.3f\n', ...
            100*it/nOut, Ebud(it,1), Ebud(it,2), Ebud(it,3));
    end
end
toc

%% 7. Checks
fprintf('\n--- checks ---\n');
fprintf('alive: %d/%d\n', nnz(P.alive), N);
nSinWater = nnz(P.matIdx==1 & ~P.p);
fprintf('S particles in the water (must be 0): %d\n', nSinWater);
fprintf('particles above the free surface (must be 0): %d\n', nnz(P.x(:,3) < -1e-9));
fprintf('energy budget closes: %.10f\n', sum(Ebud(end,:)));
fprintf('final split: water %.3f | composite P %.3f | composite S %.3f\n', Ebud(end,:));
fprintf('equipartition (vp/vs)^3 = %.3f, measured S/P in composite = %.3f\n', ...
    (comp.vp/comp.vs)^3, Ebud(end,3)/max(Ebud(end,2),eps));

%% 8. Energy budget vs time
figure('Name','Energy budget','Color','w','Position',[80 80 780 480]);
plot(1e6*tVec, Ebud(:,1), 'b-', 'LineWidth',1.6); hold on
plot(1e6*tVec, Ebud(:,2), 'r-', 'LineWidth',1.6);
plot(1e6*tVec, Ebud(:,3), 'g-', 'LineWidth',1.6);
plot(1e6*tVec, sum(Ebud,2), 'k--', 'LineWidth',1.1);
xline(1e6*(H-zsrc)/water.vp, 'k:', 'first arrival at interface', 'LineWidth',1.2);
xlabel('t (\mus)'); ylabel('energy fraction'); grid on; box on
legend('water (P only)','composite P','composite S','total','Location','best');
title('Where the source energy goes');

%% 9. Depth profile of the composite field
figure('Name','Depth profiles','Color','w','Position',[80 80 900 380]);
snap = unique(max(1, round([0.25 0.5 0.75 1]*nOut)));
lbl = arrayfun(@(k) sprintf('%.0f \\mus', 1e6*tVec(k)), snap, 'uni', 0);
nmPS = {'P','S'};
for k = 1:2
    subplot(1,2,k);
    plot(1e3*zc, squeeze(Ez(snap,:,k)).', 'LineWidth',1.2);
    xline(1e3*H, 'k--', 'water / composite', 'LineWidth',1.2);
    xlabel('depth (mm)'); ylabel('energy fraction per bin'); grid on
    title(nmPS{k});
    if k==1, legend(lbl,'Location','best'); end
end

%% 10. Surf: P and S energy vs depth and time, side by side
% The white line is the water/composite boundary. Above it the S panel is
% empty by construction -- a fluid carries no shear -- so the line is also a
% running check that nothing leaks across as the wrong mode.
figure('Name','P/S energy vs depth and time','Color','w','Position',[80 80 1220 500]);
for k = 1:2
    subplot(1,2,k);
    surf(1e3*zc, 1e6*tVec, Ez(:,:,k), 'EdgeColor','none');
    view(2); axis tight
    xlabel('depth (mm)'); ylabel('t (\mus)'); title(nmPS{k});
    colormap(gca,turbo); colorbar
    hold on
    zTop = max(Ez(:,:,k),[],'all') + eps;
    plot3(1e3*H*[1 1], 1e6*[tVec(1) tVec(end)], [zTop zTop], 'w-', 'LineWidth',2);
    hold off
end
sgtitle('white line = water / composite interface');

%% 11. Receiver array: backscattering from the heterogeneous medium
% 128 sensors on the surface, 64 either side of the source. They are in the
% water, so they record P only. The direct arrival and the specular bounce
% off the interface come first; everything after is energy the composite
% scattered back up, which is what this array is for.
tDirect = zsrc/water.vp;                        % source to surface
tSpec   = (2*H - zsrc)/water.vp;                % surface after one interface bounce

figure('Name','Backscatter B-scan','Color','w','Position',[80 80 1220 520]);

subplot(1,2,1);
Eplot = Esens;  Eplot(Eplot<=0) = NaN;          % log scale: blanks the empty cells
imagesc(1e3*xs, 1e6*tVec, log10(Eplot));
set(gca,'YDir','normal'); axis tight
xlabel('sensor position x (mm)'); ylabel('t (\mus)');
title('log_{10} energy per sensor');
colormap(gca,turbo); cb = colorbar; cb.Label.String = 'log_{10}(E/E_0)';
hold on
yline(1e6*tSpec, 'w--', 'specular bounce off the interface', 'LineWidth',1.4);
hold off

subplot(1,2,2);
semilogy(1e6*tVec, sum(Esens,2), 'k-', 'LineWidth',1.5); hold on
semilogy(1e6*tVec, sum(Esens(:,abs(xs)<=0.02),2), 'r-', 'LineWidth',1.2);
semilogy(1e6*tVec, sum(Esens(:,abs(xs)> 0.04),2), 'b-', 'LineWidth',1.2);
xline(1e6*tSpec, 'k--', 'LineWidth',1.2);
xlabel('t (\mus)'); ylabel('energy on the array'); grid on
legend('all 128 sensors','|x| <= 20 mm (near offset)','|x| > 40 mm (far offset)', ...
       'Location','best');
title('backscattered energy vs time');

fprintf('\n--- receiver array ---\n');
fprintf('128 sensors, pitch %.1f mm, span %.0f to %.0f mm, cells 0 <= z <= %.0f mm, |y| < %.0f mm\n', ...
    1e3*dxs, 1e3*xs(1), 1e3*xs(end), 1e3*zsens, 1e3*ysl);
fprintf('direct arrival at %.2f us, specular interface bounce at %.2f us\n', ...
    1e6*tDirect, 1e6*tSpec);
lateWin = tVec > 2*tSpec;
fprintf('energy on the array after 2x the specular time: %.4e of the source\n', ...
    mean(sum(Esens(lateWin,:),2)));
fprintf('  (that late coda is the backscatter the sphere medium produced)\n');

%% 12. Animation of the x-z slice
figure('Name','P/S energy, x-z slice','Color','w','Position',[80 80 1100 520]);
cl = [0 max(frames(:))];  if cl(2)<=0, cl(2)=1; end
nm = {'P','S'};
cl = [1e-8 1e-3];
for it = 1:nOut
    for k = 1:2
        subplot(1,2,k);
        imagesc(1e3*xg, 1e3*zg, squeeze(frames(it,:,:,k)).');
        set(gca,'YDir','reverse'); axis equal tight; clim(cl);
        colormap(gca,turbo); colorbar
        set(gca,"ColorScale",'log')
        hold on
        yline(0,     'w-',  'LineWidth',1.5);
        yline(1e3*H, 'w--', 'LineWidth',1.5);
        hold off
        xlabel('x (mm)'); ylabel('depth (mm)'); title(nm{k});
    end
    sgtitle(sprintf('t = %.1f \\mus   (frame %d/%d)   solid = free surface, dashed = water/composite', ...
        1e6*tVec(it), it, nOut));
    drawnow
end

%%
save("WaterLayerOverAmonodisperse.mat")
load("WaterLayerOverAmonodisperse.mat")