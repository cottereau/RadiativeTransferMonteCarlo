function allPass = FluidSolid_Pincidence()
%FLUIDSOLID_PINCIDENCE  P wave from the FLUID side hitting a solid weld.
%
%   The reverse direction of SolidFluid_Pincidence.m. ZoeppritzFluid.m is
%   only set up for solid-incident waves (see MaterialClass.ZoeppritzFluidIncident's
%   header for why the fluid-incident matrix isn't just a relabelling);
%   ZoeppritzFluidIncident builds the fluid-incident coefficients from
%   ZoeppritzFluid's validated solid-incident output via energy-flux
%   reciprocity instead. Two checks here: (1) the usual MC-vs-table
%   comparison, (2) at normal incidence, the reflected energy is checked
%   against the elementary acoustic-impedance formula
%   ((Z_solid-Z_fluid)/(Z_solid+Z_fluid))^2 -- a closed form independent of
%   both ZoeppritzFluid and the reciprocity construction, so a real
%   external check rather than only internal consistency.
%
%   See also MATERIALINTERFACESOLIDFLUID, MATERIALCLASS.ZOEPPRITZFLUIDINCIDENT.

fprintf('\n=== FluidSolid_Pincidence ===\n');
allPass = true;

vpS=3000; vsS=1500; rhoS=2200; vpF=1500; rhoF=1000;
N = 200000;
j1_deg = linspace(0,90,361);

matF = MaterialClass; matF.acoustics=false; matF.vp=vpF; matF.vs=0; matF.rho=rhoF;
zoepFluid = MaterialClass.ZoeppritzFluidIncident(j1_deg, vpS,vsS,rhoS, vpF,rhoF);

mcAng = {[],[],[]}; mcVal = {[],[],[]};
for theta = [0 15 30 45 60]
    P.x = zeros(N,3);
    P.dir = repmat([sind(theta) 0 cosd(theta)],N,1);
    P.perp = repmat([cosd(theta) 0 -sind(theta)],N,1);
    P.p = true(N,1); P.N = N;
    n = repmat([0 0 1],N,1);
    o = matF.MaterialInterfaceSolidFluid(true, zoepFluid, P, n, true(N,1));

    R   = interp1(j1_deg,zoepFluid.E_R,  theta);
    Tp  = interp1(j1_deg,zoepFluid.E_Tp, theta);
    Tsv = interp1(j1_deg,zoepFluid.E_Tsv,theta);

    fprintf('  theta = %d deg (sum theory = %.4f):\n', theta, R+Tp+Tsv);
    tol = 5/sqrt(N);
    allPass = checkChannel('R (fluid)', R,   mean(o.p2p),   tol) && allPass;
    allPass = checkChannel('T_p',       Tp,  mean(o.p2Tp),  tol) && allPass;
    allPass = checkChannel('T_sv',      Tsv, mean(o.p2Tsv), tol) && allPass;

    mcAng{1}(end+1)=theta; mcVal{1}(end+1)=mean(o.p2p);
    mcAng{2}(end+1)=theta; mcVal{2}(end+1)=mean(o.p2Tp);
    mcAng{3}(end+1)=theta; mcVal{3}(end+1)=mean(o.p2Tsv);
end

plotChannels('FluidSolid\_Pincidence', j1_deg, ...
    {zoepFluid.E_R, zoepFluid.E_Tp, zoepFluid.E_Tsv}, ...
    {'R (fluid)','T_p','T_{sv}'}, mcAng, mcVal);

% --- RTE phase: the fluid is medium 1 here, so the packet starts in it ---
rteMats = {makeInterfaceMaterial(vpF,0,rhoF), makeInterfaceMaterial(vpS,vsS,rhoS)};
expFun = @(th) [interp1(j1_deg,zoepFluid.E_R,th), NaN, NaN, ...
                interp1(j1_deg,zoepFluid.E_Tp,th), ...
                interp1(j1_deg,zoepFluid.E_Tsv,th), NaN, 1];
allPass = runEnergyPhase('FluidSolid_Pincidence', rteMats, 'P', ...
    [15 30 45 60], expFun, [1 4 5]) && allPass;

% independent check: normal-incidence acoustic-impedance formula
Zs = rhoS*vpS; Zf = rhoF*vpF;
R_imp = ((Zs-Zf)/(Zs+Zf))^2;
R_zoep = interp1(j1_deg,zoepFluid.E_R,0);
fprintf('  normal incidence: R (reciprocity) = %.6f, R (impedance formula) = %.6f\n', R_zoep, R_imp);
ok = abs(R_zoep-R_imp) < 1e-6;
allPass = ok && allPass;
fprintf('  [%s]\n', mat2char(ok));

if allPass, fprintf('FluidSolid_Pincidence: ALL PASS\n'); else, fprintf('FluidSolid_Pincidence: SOME FAILED\n'); end
end

function s = mat2char(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
