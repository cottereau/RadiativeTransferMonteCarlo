function allPass = TwoSolid_Pincidence()
%TWOSOLID_PINCIDENCE  P wave hitting a solid-solid weld: MC vs Zoeppritz.
%
%   Shoots a P wave at a flat interface between two homogeneous solids and
%   checks that the fraction of particles MaterialInterfaceTwoSolid sends
%   into each of the four channels (reflected P, reflected SV, transmitted
%   P, transmitted SV) matches ZoeppritzSolid's energy coefficients, at
%   several incidence angles for two material pairs -- one with no
%   critical angle in range, one that crosses the P->P critical angle.
%
%   See also MATERIALINTERFACETWOSOLID, ZOEPPRITZSOLID.

fprintf('\n=== TwoSolid_Pincidence ===\n');
allPass = true;

pairs = {
    struct('name','shale/sandstone (no critical angle in range)', ...
           'vp1',3000,'vs1',1500,'rho1',2200,'vp2',4500,'vs2',2600,'rho2',2600, ...
           'theta',[10 25 45 60])
    struct('name','slow/fast (crosses the P->P critical angle)', ...
           'vp1',1600,'vs1',800,'rho1',1800,'vp2',3500,'vs2',2000,'rho2',2600, ...
           'theta',[10 20 35 50])
};

N = 200000;
j1_deg = linspace(0,90,361);

for k = 1:numel(pairs)
    pr = pairs{k};
    fprintf('  -- %s --\n', pr.name);
    matA = MaterialClass; matA.acoustics=false; matA.vp=pr.vp1; matA.vs=pr.vs1; matA.rho=pr.rho1;
    matB = MaterialClass; matB.acoustics=false; matB.vp=pr.vp2; matB.vs=pr.vs2; matB.rho=pr.rho2;
    zoep = MaterialClass.ZoeppritzSolid(j1_deg,matA.vp,matA.vs,matA.rho,matB.vp,matB.vs,matB.rho);
    zoep.j1_deg = j1_deg;

    mcAng = {[],[],[],[]}; mcVal = {[],[],[],[]};
    for theta = pr.theta
        P.x = zeros(N,3);
        P.dir = repmat([sind(theta) 0 cosd(theta)],N,1);
        P.perp = repmat([cosd(theta) 0 -sind(theta)],N,1);
        P.p = true(N,1); P.N = N;
        n = repmat([0 0 1],N,1);
        o = matA.MaterialInterfaceTwoSolid(matB, zoep, P, n, true(N,1));

        Rpp  = interp1(j1_deg,zoep.E_Rpp, theta);
        Rpsv = interp1(j1_deg,zoep.E_Rpsv,theta);
        Tpp  = interp1(j1_deg,zoep.E_Tpp, theta);
        Tpsv = interp1(j1_deg,zoep.E_Tpsv,theta);

        fprintf('  theta = %d deg (sum theory = %.4f):\n', theta, Rpp+Rpsv+Tpp+Tpsv);
        tol = 5/sqrt(N); % ~5 sigma of a binomial-fraction estimate at N particles
        allPass = checkChannel('R_pp',  Rpp,  mean(o.p2p),  tol) && allPass;
        allPass = checkChannel('R_psv', Rpsv, mean(o.p2sv), tol) && allPass;
        allPass = checkChannel('T_pp',  Tpp,  mean(o.p2Tp), tol) && allPass;
        allPass = checkChannel('T_psv', Tpsv, mean(o.p2Tsv),tol) && allPass;

        mcAng{1}(end+1)=theta; mcVal{1}(end+1)=mean(o.p2p);
        mcAng{2}(end+1)=theta; mcVal{2}(end+1)=mean(o.p2sv);
        mcAng{3}(end+1)=theta; mcVal{3}(end+1)=mean(o.p2Tp);
        mcAng{4}(end+1)=theta; mcVal{4}(end+1)=mean(o.p2Tsv);
    end

    plotChannels(sprintf('TwoSolid\\_Pincidence: %s', pr.name), j1_deg, ...
        {zoep.E_Rpp, zoep.E_Rpsv, zoep.E_Tpp, zoep.E_Tpsv}, ...
        {'R_{pp}','R_{psv}','T_{pp}','T_{psv}'}, mcAng, mcVal);

    % --- RTE phase: same pair and angles, but the transport actually runs ---
    rteMats = {makeInterfaceMaterial(pr.vp1,pr.vs1,pr.rho1), ...
               makeInterfaceMaterial(pr.vp2,pr.vs2,pr.rho2)};
    expFun = @(th) [interp1(j1_deg,zoep.E_Rpp,th), interp1(j1_deg,zoep.E_Rpsv,th), NaN, ...
                    interp1(j1_deg,zoep.E_Tpp,th), interp1(j1_deg,zoep.E_Tpsv,th), NaN, 1];
    allPass = runEnergyPhase(sprintf('TwoSolid_Pincidence (%s)',pr.name), ...
        rteMats, 'P', pr.theta, expFun, [1 2 4 5]) && allPass;
end

if allPass, fprintf('TwoSolid_Pincidence: ALL PASS\n'); else, fprintf('TwoSolid_Pincidence: SOME FAILED\n'); end
end
