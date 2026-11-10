function deriveSphereBC(outFile)
%% deriveSphereBC
% Symbolically derives the six boundary quantities
% (u_r, u_th, u_ph, sigma_rr, sigma_rth, sigma_rph) of the three basic
% partial-wave fields of a sphere, for any azimuthal order m, and generates
% the numerical kernel used by TMatrix.sphereScattering.
%
%   TMatrix.deriveSphereBC()    % writes sphereBCkernel.m next to this file
%
% Method: everything in Cartesian coordinates, where grad, curl and div have
% no curvilinear-coordinate pitfalls; only at the end is the result projected
% on r_hat, th_hat and ph_hat. The only things written "from memory" are the
% three potential definitions:
%
%   P esferoidal  :  u = grad(phi)              phi = Z(r) L(cos th) cos(m ph)
%   SV esferoidal :  u = curl(curl(r*chi))      chi = Z(r) L(cos th) cos(m ph)
%   SH toroidal   :  u = curl(r*psi)            psi = Z(r) L(cos th) sin(m ph)
%
% Z and L enter as cubic Taylor expansions around the evaluation point (exact
% up to the third derivative, the highest order the stress needs). The
% azimuthal dependence stays analytic. No field equation (Helmholtz, Legendre)
% is used, so they remain available as independent checks.
%
% m = 0 recovers the axisymmetric case: the toroidal terms and the ph
% components vanish identically.

if nargin < 1
    outFile = fullfile(fileparts(mfilename('fullpath')),'sphereBCkernel.m');
end

syms x y z real
syms lam mu real
syms R positive
syms th Ph real
syms mm real
syms Z0 Z1 Z2 Z3 L0 L1 L2 L3 real

rr = sqrt(x^2 + y^2 + z^2);
uu = z/rr;                      % cos(theta)
pp = atan2(y,x);                % phi -- atan2, not atan(y/x): with atan the
                                % kernel would only hold for |phi| < pi/2, and
                                % the violation shows up as a silent sign
                                % flip in the radial factors

dr = rr - R;
du = uu - cos(th);

Zex = Z0 + Z1*dr + Z2*dr^2/2 + Z3*dr^3/6;
Lex = L0 + L1*du + L2*du^2/2 + L3*du^3/6;

radial_ang = Zex*Lex;
Pot_c = radial_ang*cos(mm*pp);      % spheroidal
Pot_s = radial_ang*sin(mm*pp);      % toroidal

grad = @(f) [diff(f,x); diff(f,y); diff(f,z)];
divg = @(v) diff(v(1),x) + diff(v(2),y) + diff(v(3),z);
curl = @(v) [diff(v(3),y)-diff(v(2),z); ...
             diff(v(1),z)-diff(v(3),x); ...
             diff(v(2),x)-diff(v(1),y)];

rvec = [x; y; z];

uP = grad(Pot_c);                       % P spheroidal
uS = curl(curl(rvec*Pot_c));            % SV spheroidal
uT = curl(rvec*Pot_s);                  % SH toroidal

tensao = @(u) lam*divg(u)*eye(3) + mu*(jacobian(u,[x y z]) + jacobian(u,[x y z]).');

rhat  = [sin(th)*cos(Ph); sin(th)*sin(Ph); cos(th)];
thhat = [cos(th)*cos(Ph); cos(th)*sin(Ph); -sin(th)];
phhat = [-sin(Ph); cos(Ph); 0];
sub0  = @(e) subs(e,[x y z],(R*rhat).');

bc = sym(zeros(18,1));
campos = {uP,uS,uT};
for ic = 1:3
    u = campos{ic};
    S = tensao(u);
    k0 = 6*(ic-1);
    bc(k0+1) = sub0(u.'*rhat);
    bc(k0+2) = sub0(u.'*thhat);
    bc(k0+3) = sub0(u.'*phhat);
    bc(k0+4) = sub0(rhat.'*S*rhat);
    bc(k0+5) = sub0(thhat.'*S*rhat);
    bc(k0+6) = sub0(phhat.'*S*rhat);
end
bc = simplify(bc,'Steps',20);

allowed = [R th Ph lam mu mm Z0 Z1 Z2 Z3 L0 L1 L2 L3];
leftover = setdiff(symvar(bc),allowed);
if ~isempty(leftover)
    error('unexpected leftover symbol: %s',char(leftover(1)));
end

matlabFunction(bc,'File',outFile, ...
    'Vars',{R,th,Ph,mm,lam,mu,Z0,Z1,Z2,Z3,L0,L1,L2,L3}, ...
    'Outputs',{'bc'});

fprintf('generated %s\n',outFile);
fprintf(['bc(1:6)   = P spheroidal field  [ur uth uph srr srt srp]\n' ...
         'bc(7:12)  = SV spheroidal field\n' ...
         'bc(13:18) = SH toroidal field\n']);
end
