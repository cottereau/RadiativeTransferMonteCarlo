function allPass = SolidFluid_Pincidence()
%SOLIDFLUID_PINCIDENCE  P wave from the SOLID side hitting a fluid weld.
%
%   A fluid has no S channel, so this is a three-way check: reflected P,
%   reflected SV (mode-converted, stays in the solid) and transmitted P
%   (into the fluid -- there is no transmitted SV, a fluid cannot carry
%   one) against ZoeppritzFluid's P-incidence energy coefficients, for a
%   water-like and a gas-like fluid.
%
%   See also MATERIALINTERFACESOLIDFLUID, ZOEPPRITZFLUID.

fprintf('\n=== SolidFluid_Pincidence ===\n');
allPass = true;

fluids = {
    struct('name','sandstone / water', 'vpF',1500,'rhoF',1000,'theta',[10 25 40 60])
    struct('name','sandstone / gas (near-vacuum limit)', 'vpF',340,'rhoF',1.2,'theta',[0 15 30])
};
vpS=3000; vsS=1500; rhoS=2200;
N = 200000;
j1_deg = linspace(0,90,361);

for k = 1:numel(fluids)
    fl = fluids{k};
    fprintf('  -- %s --\n', fl.name);

    matS = MaterialClass; 
    matS.acoustics=false; 
    matS.vp=vpS; 
    matS.vs=vsS; 
    matS.rho=rhoS;

    matF = MaterialClass; 
    matF.acoustics=false; 
    matF.vp=fl.vpF; 
    matF.vs=0; 
    matF.rho=fl.rhoF;

    zoepSolid.P  = MaterialClass.ZoeppritzFluid(j1_deg, vpS,vsS,rhoS, fl.vpF,0,fl.rhoF);
    zoepSolid.SV = zoepSolid.P;
    zoepSolid.j1_deg = j1_deg;
    zoepSolid.vpFluid = fl.vpF;

    mcAng = {[],[],[]}; mcVal = {[],[],[]};
    for theta = fl.theta
        P.x = zeros(N,3);
        P.dir = repmat([sind(theta) 0 cosd(theta)],N,1);
        P.perp = repmat([cosd(theta) 0 -sind(theta)],N,1);
        P.p = true(N,1); P.N = N;
        n = repmat([0 0 1],N,1);
        o = matS.MaterialInterfaceSolidFluid(false, zoepSolid, P, n, true(N,1));

        Rpp  = interp1(j1_deg,zoepSolid.P.E_Rpp, theta);
        Rpsv = interp1(j1_deg,zoepSolid.P.E_Rpsv,theta);
        Tpp  = interp1(j1_deg,zoepSolid.P.E_Tpp, theta);

        fprintf('  theta = %d deg (sum theory = %.4f):\n', theta, Rpp+Rpsv+Tpp);
        tol = 5/sqrt(N);
        allPass = checkChannel('R_pp',  Rpp,  mean(o.p2p),  tol) && allPass;
        allPass = checkChannel('R_psv', Rpsv, mean(o.p2sv), tol) && allPass;
        allPass = checkChannel('T_p',   Tpp,  mean(o.p2Tp), tol) && allPass;

        mcAng{1}(end+1)=theta; mcVal{1}(end+1)=mean(o.p2p);
        mcAng{2}(end+1)=theta; mcVal{2}(end+1)=mean(o.p2sv);
        mcAng{3}(end+1)=theta; mcVal{3}(end+1)=mean(o.p2Tp);
    end

    plotChannels(sprintf('SolidFluid\\_Pincidence: %s', fl.name), j1_deg, ...
        {zoepSolid.P.E_Rpp, zoepSolid.P.E_Rpsv, zoepSolid.P.E_Tpp}, ...
        {'R_{pp}','R_{psv}','T_p'}, mcAng, mcVal);

    % --- RTE phase: same pair and angles, but the transport actually runs ---
    rteMats = {makeInterfaceMaterial(vpS,vsS,rhoS), ...
               makeInterfaceMaterial(fl.vpF,0,fl.rhoF)};
    expFun = @(th) [interp1(j1_deg,zoepSolid.P.E_Rpp,th), ...
                    interp1(j1_deg,zoepSolid.P.E_Rpsv,th), NaN, ...
                    interp1(j1_deg,zoepSolid.P.E_Tpp,th), NaN, NaN, 1];
    allPass = runEnergyPhase(sprintf('SolidFluid_Pincidence (%s)',fl.name), ...
        rteMats, 'P', fl.theta, expFun, [1 2 4]) && allPass;
end

if allPass, fprintf('SolidFluid_Pincidence: ALL PASS\n'); else, fprintf('SolidFluid_Pincidence: SOME FAILED\n'); end
end
