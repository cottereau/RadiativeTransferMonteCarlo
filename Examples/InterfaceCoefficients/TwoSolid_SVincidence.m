function allPass = TwoSolid_SVincidence()
%TWOSOLID_SVINCIDENCE  SV wave hitting a solid-solid weld: MC vs Zoeppritz.
%
%   Same idea as TwoSolid_Pincidence.m, but for an SV wave: checks
%   reflected SV, reflected P (mode-converted), transmitted SV and
%   transmitted P against ZoeppritzSolid's SV-incidence energy
%   coefficients. The SV polarization is set in the plane of incidence
%   (P.perp coplanar with dir and the normal) so MaterialInterfaceTwoSolid
%   classifies every particle as SV, not SH -- see
%   TwoSolid_SHincidence.m for the orthogonal case.
%
%   See also MATERIALINTERFACETWOSOLID, ZOEPPRITZSOLID.

fprintf('\n=== TwoSolid_SVincidence ===\n');
allPass = true;

pairs = {
    struct('name','shale/sandstone', ...
           'vp1',3000,'vs1',1500,'rho1',2200,'vp2',4500,'vs2',2600,'rho2',2600, ...
           'theta',[10 20 30 40])
    struct('name','slow/fast', ...
           'vp1',1600,'vs1',800,'rho1',1800,'vp2',3500,'vs2',2000,'rho2',2600, ...
           'theta',[10 20 30 40])
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
        P.perp = repmat([cosd(theta) 0 -sind(theta)],N,1); % in-plane: pure SV
        P.p = false(N,1); P.N = N;
        n = repmat([0 0 1],N,1);
        o = matA.MaterialInterfaceTwoSolid(matB, zoep, P, n, true(N,1));

        Rsvsv = interp1(j1_deg,zoep.E_Rsvsv,theta);
        Rsp   = interp1(j1_deg,zoep.E_Rsp,  theta);
        Tsvsv = interp1(j1_deg,zoep.E_Tsvsv,theta);
        Tsp   = interp1(j1_deg,zoep.E_Tsp,  theta);

        fprintf('  theta = %d deg (sum theory = %.4f):\n', theta, Rsvsv+Rsp+Tsvsv+Tsp);
        tol = 5/sqrt(N);
        allPass = checkChannel('R_svsv', Rsvsv, mean(o.sv2sv), tol) && allPass;
        allPass = checkChannel('R_sp',   Rsp,   mean(o.sv2p),  tol) && allPass;
        allPass = checkChannel('T_svsv', Tsvsv, mean(o.sv2Tsv),tol) && allPass;
        allPass = checkChannel('T_sp',   Tsp,   mean(o.sv2Tp), tol) && allPass;

        mcAng{1}(end+1)=theta; mcVal{1}(end+1)=mean(o.sv2sv);
        mcAng{2}(end+1)=theta; mcVal{2}(end+1)=mean(o.sv2p);
        mcAng{3}(end+1)=theta; mcVal{3}(end+1)=mean(o.sv2Tsv);
        mcAng{4}(end+1)=theta; mcVal{4}(end+1)=mean(o.sv2Tp);
    end

    plotChannels(sprintf('TwoSolid\\_SVincidence: %s', pr.name), j1_deg, ...
        {zoep.E_Rsvsv, zoep.E_Rsp, zoep.E_Tsvsv, zoep.E_Tsp}, ...
        {'R_{svsv}','R_{sp}','T_{svsv}','T_{sp}'}, mcAng, mcVal);

    % --- RTE phase: same pair and angles, but the transport actually runs ---
    rteMats = {makeInterfaceMaterial(pr.vp1,pr.vs1,pr.rho1), ...
               makeInterfaceMaterial(pr.vp2,pr.vs2,pr.rho2)};
    expFun = @(th) [interp1(j1_deg,zoep.E_Rsp,th), interp1(j1_deg,zoep.E_Rsvsv,th), NaN, ...
                    interp1(j1_deg,zoep.E_Tsp,th), interp1(j1_deg,zoep.E_Tsvsv,th), NaN, 1];
    allPass = runEnergyPhase(sprintf('TwoSolid_SVincidence (%s)',pr.name), ...
        rteMats, 'SV', pr.theta, expFun, [1 2 4 5]) && allPass;
end

if allPass, fprintf('TwoSolid_SVincidence: ALL PASS\n'); else, fprintf('TwoSolid_SVincidence: SOME FAILED\n'); end
end
