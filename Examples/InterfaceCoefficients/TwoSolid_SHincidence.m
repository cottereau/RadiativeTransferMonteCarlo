function allPass = TwoSolid_SHincidence()
%TWOSOLID_SHINCIDENCE  SH wave hitting a solid-solid weld: MC vs Zoeppritz.
%
%   SH never mode-converts to P/SV at a solid-solid interface (it stays SH
%   on both sides), so this is a simpler two-channel check: reflected SH
%   and transmitted SH against ZoeppritzSolid's E_Rsh/E_Tsh. The
%   polarization is set orthogonal to the plane of incidence (out of the
%   x-z plane here) so MaterialInterfaceTwoSolid classifies every particle
%   as SH.
%
%   See also MATERIALINTERFACETWOSOLID, ZOEPPRITZSOLID.

fprintf('\n=== TwoSolid_SHincidence ===\n');
allPass = true;

pr = struct('vp1',3000,'vs1',1500,'rho1',2200,'vp2',4500,'vs2',2600,'rho2',2600);
N = 200000;
j1_deg = linspace(0,90,361);

matA = MaterialClass; matA.acoustics=false; matA.vp=pr.vp1; matA.vs=pr.vs1; matA.rho=pr.rho1;
matB = MaterialClass; matB.acoustics=false; matB.vp=pr.vp2; matB.vs=pr.vs2; matB.rho=pr.rho2;
zoep = MaterialClass.ZoeppritzSolid(j1_deg,matA.vp,matA.vs,matA.rho,matB.vp,matB.vs,matB.rho);
zoep.j1_deg = j1_deg;

mcAng = {[],[]}; mcVal = {[],[]};
for theta = [10 20 35 50 65]
    P.x = zeros(N,3);
    P.dir = repmat([sind(theta) 0 cosd(theta)],N,1);
    P.perp = repmat([0 1 0],N,1); % out of the plane of incidence: pure SH
    P.p = false(N,1); P.N = N;
    n = repmat([0 0 1],N,1);
    o = matA.MaterialInterfaceTwoSolid(matB, zoep, P, n, true(N,1));

    Rsh = interp1(j1_deg,zoep.E_Rsh,theta);
    Tsh = interp1(j1_deg,zoep.E_Tsh,theta);

    fprintf('  theta = %d deg (sum theory = %.4f, frac classified SH = %.4f):\n', ...
        theta, Rsh+Tsh, nnz(o.shparticle)/nnz(o.sparticle));
    tol = 5/sqrt(N);
    allPass = checkChannel('R_sh', Rsh, mean(o.sh2sh),  tol) && allPass;
    allPass = checkChannel('T_sh', Tsh, mean(o.sh2Tsh), tol) && allPass;

    mcAng{1}(end+1)=theta; mcVal{1}(end+1)=mean(o.sh2sh);
    mcAng{2}(end+1)=theta; mcVal{2}(end+1)=mean(o.sh2Tsh);
end

plotChannels('TwoSolid\_SHincidence', j1_deg, {zoep.E_Rsh, zoep.E_Tsh}, ...
    {'R_{sh}','T_{sh}'}, mcAng, mcVal);

% --- RTE phase: same pair and angles, but the transport actually runs ---
rteMats = {makeInterfaceMaterial(pr.vp1,pr.vs1,pr.rho1), ...
           makeInterfaceMaterial(pr.vp2,pr.vs2,pr.rho2)};
expFun = @(th) [NaN, NaN, interp1(j1_deg,zoep.E_Rsh,th), ...
                NaN, NaN, interp1(j1_deg,zoep.E_Tsh,th), 1];
allPass = runEnergyPhase('TwoSolid_SHincidence', rteMats, 'SH', ...
    [10 20 35 50 65], expFun, [3 6]) && allPass;

if allPass, fprintf('TwoSolid_SHincidence: ALL PASS\n'); else, fprintf('TwoSolid_SHincidence: SOME FAILED\n'); end
end
