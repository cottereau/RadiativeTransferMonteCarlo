function obs = HomogeneoCorrigido(numParticles,outFile)
%HOMOGENEOCORRIGIDO  Homogeneo.m with the corrected initial condition.
%
%   Same cylinder, material and observation grid as Homogeneo.m, with three
%   corrections:
%     1. the radial profile is rebuilt from the SEM source time function as
%        Mddot(t0-r/alpha)^2 / r^2 (see initialConditionRTE), instead of the
%        misplaced ICrteManda2 shell;
%     2. directions are 'outgoing': a free wave packet from a point source
%        is purely radial;
%     3. source.halfspace folds the particles born above the free surface
%        back into the solid rather than losing them, which with the
%        corrected shell (radius ~20 mm) would throw away half the energy.
%
%   The RTE clock then starts at t0, so traces compare to SEM time minus t0.

if nargin < 1 || isempty(numParticles), numParticles = 5e6; end
if nargin < 2 || isempty(outFile), outFile = 'HomogeneoCorrigido.mat'; end

base  = '/Users/lac/Nextcloud/ArtigoManda';
t0    = 6e-6;
vvp   = 4814.4;  vvs = 2962.1;  rho = 2210.1;  freq = 500e3;

IC = initialConditionRTE([base '/SEM/CasoHomogeneoEInicial/Homo/AmpTab.txt'],vvp,t0);
fprintf('CI: casca de %.2f a %.2f mm, pico em %.2f mm (t0 = %.1f us)\n', ...
    1e3*IC.rInner,1e3*IC.rOuter,1e3*IC.rPeak,1e6*t0);

geometry = struct('dimension',3,'frame','cylindrical');
geometry.bnd(1) = struct('dir',3,'val',0,'axis',0);
geometry.bnd(2) = struct('dir',3,'val',-0.23,'axis',0);
geometry.bnd(3) = struct('dir',4,'val',0.055,'axis',3);

source = struct('numberParticles',numParticles, ...
                'type','point', ...
                'position',[0 0 0], ...
                'direction','outgoing', ...
                'radial',griddedInterpolant(IC.r,IC.amp), ...
                'halfspace',[3 0 +1], ...      % fold z > 0 back into z < 0
                'polarization','P');           % isotropic moment radiates P only

xx = [0 2.5e-3:5e-4:5.5e-2];
observation = struct('x',xx,'y',[-pi pi],'z',-0.23:5e-4:0, ...
                     'directions',[0 pi],'time',0:0.5e-6:2e-4);

mat = MaterialClass;
mat.acoustics = false;
mat.vp = vvp;  mat.vs = vvs;  mat.rho = rho;  mat.Frequency = freq;
% Homogeneous medium. Homogeneo.m sets the coefficients of variation to
% 1e-100 and lets CalcSigma build the cross sections from the PSDF, which
% underflows to a degenerate CDF in the current MaterialClass. Supplying
% negligible isotropic cross sections directly is the same physics and is
% numerically well behaved: prepareSigma skips CalcSigma when sigma is set.
tiny = @(th) 1e-30*ones(size(th));
mat.sigma = {tiny tiny; tiny tiny};

material = mat.copyobj;
material.timeSteps = 0;
material = MaterialClass.prepareSigma(material,3);
fprintf('livre caminho medio: %s m (corpo tem 0.23 m)\n', ...
    mat2str(material.meanFreePath(:).',3));

tic; obs = radiativeTransfer(geometry,source,material,observation); toc

dV  = (obs.x(:).*diff(obs.binX(:)))*(diff(obs.binZ(:)).')*diff(obs.binY);
tot = squeeze(sum(sum(obs.energyDensity.*dV,1),2));
fprintf('energia dentro do corpo: %.4f (era 0.4999 sem o halfspace)\n', ...
    median(sum(tot(obs.t > 2e-5,:),2)));

save(outFile,'obs','IC','t0','-v7.3');
fprintf('salvo em %s\n',outFile);
end
