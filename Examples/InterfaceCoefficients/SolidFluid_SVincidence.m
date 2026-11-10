function allPass = SolidFluid_SVincidence()
%SOLIDFLUID_SVINCIDENCE  SV wave from the SOLID side hitting a fluid weld.
%
%   Reflected SV, reflected P (mode-converted) and transmitted P (into the
%   fluid, mode-converted -- a fluid only ever carries P) against
%   ZoeppritzFluid's SV-incidence energy coefficients.
%
%   See also MATERIALINTERFACESOLIDFLUID, ZOEPPRITZFLUID.

fprintf('\n=== SolidFluid_SVincidence ===\n');
allPass = true;

vpS=3000; vsS=1500; rhoS=2200; vpF=1500; rhoF=1000;
N = 200000;
j1_deg = linspace(0,90,361);

matS = MaterialClass; matS.acoustics=false; matS.vp=vpS; matS.vs=vsS; matS.rho=rhoS;
zoepSolid.P  = MaterialClass.ZoeppritzFluid(j1_deg, vpS,vsS,rhoS, vpF,0,rhoF);
zoepSolid.SV = zoepSolid.P;
zoepSolid.j1_deg = j1_deg;
zoepSolid.vpFluid = vpF;

mcAng = {[],[],[]}; mcVal = {[],[],[]};
for theta = [10 20 30 40]
    P.x = zeros(N,3);
    P.dir = repmat([sind(theta) 0 cosd(theta)],N,1);
    P.perp = repmat([cosd(theta) 0 -sind(theta)],N,1); % in-plane: pure SV
    P.p = false(N,1); P.N = N;
    n = repmat([0 0 1],N,1);
    o = matS.MaterialInterfaceSolidFluid(false, zoepSolid, P, n, true(N,1));

    Rsvsv = interp1(j1_deg,zoepSolid.SV.E_Rsvsv,theta);
    Rsp   = interp1(j1_deg,zoepSolid.SV.E_Rsp,  theta);
    Tsp   = interp1(j1_deg,zoepSolid.SV.E_Tsp,  theta);

    fprintf('  theta = %d deg (sum theory = %.4f):\n', theta, Rsvsv+Rsp+Tsp);
    tol = 5/sqrt(N);
    allPass = checkChannel('R_svsv', Rsvsv, mean(o.sv2sv), tol) && allPass;
    allPass = checkChannel('R_sp',   Rsp,   mean(o.sv2p),  tol) && allPass;
    allPass = checkChannel('T_p',    Tsp,   mean(o.sv2Tp), tol) && allPass;

    mcAng{1}(end+1)=theta; mcVal{1}(end+1)=mean(o.sv2sv);
    mcAng{2}(end+1)=theta; mcVal{2}(end+1)=mean(o.sv2p);
    mcAng{3}(end+1)=theta; mcVal{3}(end+1)=mean(o.sv2Tp);
end

plotChannels('SolidFluid\_SVincidence', j1_deg, ...
    {zoepSolid.SV.E_Rsvsv, zoepSolid.SV.E_Rsp, zoepSolid.SV.E_Tsp}, ...
    {'R_{svsv}','R_{sp}','T_p'}, mcAng, mcVal);

% --- RTE phase: same pair and angles, but the transport actually runs ---
rteMats = {makeInterfaceMaterial(vpS,vsS,rhoS), makeInterfaceMaterial(vpF,0,rhoF)};
expFun = @(th) [interp1(j1_deg,zoepSolid.SV.E_Rsp,th), ...
                interp1(j1_deg,zoepSolid.SV.E_Rsvsv,th), NaN, ...
                interp1(j1_deg,zoepSolid.SV.E_Tsp,th), NaN, NaN, 1];
allPass = runEnergyPhase('SolidFluid_SVincidence', rteMats, 'SV', ...
    [10 20 30 40], expFun, [1 2 4]) && allPass;

if allPass, fprintf('SolidFluid_SVincidence: ALL PASS\n'); else, fprintf('SolidFluid_SVincidence: SOME FAILED\n'); end
end
