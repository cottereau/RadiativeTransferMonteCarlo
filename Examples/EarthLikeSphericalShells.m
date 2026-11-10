%% Earth-like concentric spherical shells: inner core (solid) / outer core
%  (fluid) / mantle (solid), via propagateParticleMultiMaterial's dir==5
%  (spherical) interface support.
%
%  Not a PREM-style continuously-varying velocity model -- three
%  piecewise-homogeneous shells, each with its own MaterialClass, welded
%  at r=0.10 (inner/outer core, solid-fluid) and r=0.25 (outer core/mantle,
%  fluid-solid). Outer boundary r=0.40 is a reflective (free) surface.
%
%  See also PROPAGATEPARTICLEMULTIMATERIAL, MATERIALINTERFACESOLIDFLUID,
%  MATERIALINTERFACETWOSOLID.
ccc

rng(42);
matInner  = isoMat(5000,3000,3200, 0.04);   % r < 0.10          (inner core, solid)
matOuter  = isoMat(3200,0,   3000, 0.04);   % 0.10 < r < 0.25   (outer core, fluid)
matMantle = isoMat(4200,2400,2600, 0.04);   % 0.25 < r < 0.40   (mantle, solid)
matArray = {matInner, matOuter, matMantle};

geometry.dimension = 3; geometry.frame = 'cartesian';
geometry.bnd(1) = struct('dir',5,'val',0.40,'type','reflective','mat',[]); % outer free surface
geometry.bnd(2) = struct('dir',5,'val',0.10,'type','','mat',[1 2]);        % inner core / outer core
geometry.bnd(3) = struct('dir',5,'val',0.25,'type','','mat',[2 3]);        % outer core / mantle

N = 20000;
P.N = N;
P.x = zeros(N,3);                        % point source at the centre
c0 = 2*rand(N,1)-1; ph = 2*pi*rand(N,1); s0 = sqrt(1-c0.^2);
P.dir  = [s0.*cos(ph) s0.*sin(ph) c0];
P.perp = [-sin(ph) cos(ph) zeros(N,1)];
P.p = true(N,1);
P.alive = true(N,1);
P.t = zeros(N,1);
P.matIdx = ones(N,1);

T = 0.05;

% Output times. The run is stepped instead of advanced in one call so the
% P and S energy can be recorded AS it propagates -- the final state at t = T
% is the same either way, and the validation block below is unchanged.
nOut      = 120;                  % output times (animation frames)
nr        = 120;                  % radial bins for the energy history
ngrid     = 140;                  % x-z bins for the slice animation
slab      = 0.02*0.40;            % half-thickness in y of the animated slice

tVec   = linspace(0, T, nOut);
rEdges = linspace(0, 0.40, nr+1);   
rc = 0.5*(rEdges(1:end-1)+rEdges(2:end));
gEdges = linspace(-0.40, 0.40, ngrid+1); 
gc = 0.5*(gEdges(1:end-1)+gEdges(2:end));

EP     = zeros(nOut, nr);
ES     = zeros(nOut, nr);
frames = zeros(nOut, ngrid, ngrid, 2);

tic
for it = 1:nOut
    if it > 1
        P = propagateParticleMultiMaterial(matArray, geometry, P, tVec(it));
    end
    100*it/nOut
    rr  = vecnorm(P.x,2,2);
    isP = P.p & P.alive;
    isS = ~P.p & P.alive;

    EP(it,:) = histcounts(rr(isP), rEdges)/N;
    ES(it,:) = histcounts(rr(isS), rEdges)/N;

    inSlab = abs(P.x(:,2)) < slab;
    frames(it,:,:,1) = histcounts2(P.x(isP&inSlab,1), P.x(isP&inSlab,3), gEdges, gEdges)/N;
    frames(it,:,:,2) = histcounts2(P.x(isS&inSlab,1), P.x(isS&inSlab,3), gEdges, gEdges)/N;
end
toc

r = vecnorm(P.x,2,2);

fprintf('\n--- validation ---\n');
fprintf('alive: %d/%d\n', nnz(P.alive), N);
fprintf('matIdx in {1,2,3} only: %d\n', all(ismember(P.matIdx,[1 2 3])));
fprintf('fraction per shell: inner=%.4f  outer(fluid)=%.4f  mantle=%.4f\n', ...
    mean(P.matIdx==1),mean(P.matIdx==2),mean(P.matIdx==3));
fprintf('any NaN in x: %d, in dir: %d\n', any(isnan(P.x(:))), any(isnan(P.dir(:))));
fprintf('max |norm(dir)-1|: %.2e\n', max(abs(vecnorm(P.dir,2,2)-1)));
fprintf('r range: [%.4f, %.4f]  (domain [0,0.40])\n', min(r), max(r));
nBad = nnz(P.matIdx==1 & r>0.10) + nnz(P.matIdx==2 & r<0.10) + ...
       nnz(P.matIdx==2 & r>0.25) + nnz(P.matIdx==3 & r<0.25);
fprintf('cross-contamination count (should be 0): %d\n', nBad);
fprintf('outer-core (fluid) particles with S polarization (should be 0): %d\n', ...
    nnz(P.matIdx==2 & ~P.p));

figure;
edges = 0:0.005:0.40;
histogram(r,edges,'FaceColor',[0.3 0.5 0.8]);
hold on
xline(0.10,'r--','LineWidth',1.5);
xline(0.25,'r--','LineWidth',1.5);
xline(0.40,'k-','LineWidth',1.5);
xlabel('r'); ylabel('particle count');
title(sprintf('Earth-like spherical shells, radial distribution at t=%.3g', T));
legend('particles','inner/outer core','outer core/mantle','free surface','Location','best');
box on

%% P and S energy vs radius and time
% The whole history in one static figure. The fluid outer core is the point
% of interest: the S panel has to stay dark between r = 0.10 and r = 0.25,
% because that shell is fluid (vs == 0) and carries no shear wave, while the
% P panel shows energy crossing it freely.
figure('Name','P/S energy vs radius and time','Color','w','Position',[80 80 1200 470]);
modeName = {'P energy','S energy'};
modeData = {EP, ES};
for k = 1:2
    subplot(1,2,k);
    surf(rc, tVec, modeData{k}, 'EdgeColor','none');
    view(2); axis tight
    xlabel('r'); ylabel('t'); title(modeName{k});
    colormap(gca, turbo); colorbar
    hold on
    zTop = max(modeData{k}(:)) + eps;
    plot3([0.10 0.10],[tVec(1) tVec(end)],[zTop zTop],'w--','LineWidth',1.2);
    plot3([0.25 0.25],[tVec(1) tVec(end)],[zTop zTop],'w--','LineWidth',1.2);
    plot3([0.40 0.40],[tVec(1) tVec(end)],[zTop zTop],'w-','LineWidth',1.2);
    hold off
end
sgtitle('dashed = shell interfaces, solid = free surface');

fprintf('\nS energy inside the fluid shell (0.10 < r < 0.25): max = %.3e\n', ...
    max(max(ES(:, rc>0.10 & rc<0.25))));
fprintf('P energy inside the same shell:                    max = %.3e\n', ...
    max(max(EP(:, rc>0.10 & rc<0.25))));

%% Animation of the x-z slice
% Replayed on screen with drawnow. The frames are kept in `frames`, so the
% loop can be re-run at any speed, or a single instant inspected, without
% propagating again.
figAnim = figure('Name','P/S energy, x-z slice','Color','w','Position',[80 80 1100 520]);
cl = [0 max(frames(:))];
if cl(2) <= 0, cl(2) = 1; end
th = linspace(0,2*pi,200);
cl = [0 1e-3]
for it = 1:nOut
    for k = 1:2
        subplot(1,2,k);
        imagesc(gc, gc, squeeze(frames(it,:,:,k)).');
        set(gca,'YDir','normal'); axis equal tight; clim(cl);
        colormap(gca, turbo); colorbar
        xlabel('x'); ylabel('z'); title(modeName{k});
        hold on
        plot(0.10*cos(th),0.10*sin(th),'w--','LineWidth',1.0);
        plot(0.25*cos(th),0.25*sin(th),'w--','LineWidth',1.0);
        plot(0.40*cos(th),0.40*sin(th),'w-','LineWidth',1.2);
        hold off
    end
    sgtitle(sprintf('t = %.4g   (frame %d/%d)', tVec(it), it, nOut));
    drawnow
end

function mat = isoMat(vp,vs,rho,mfp)
mat = MaterialClass; mat.acoustics=false; mat.vp=vp; mat.vs=vs; mat.rho=rho;
c = 1/(4*pi*mfp);
mat.sigma = {@(th)c*ones(size(th)) @(th)c*ones(size(th)); ...
             @(th)c*ones(size(th)) @(th)c*ones(size(th))};
mat = MaterialClass.prepareSigma(mat,3);
end
