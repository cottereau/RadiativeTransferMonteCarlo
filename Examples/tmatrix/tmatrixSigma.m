function out = tmatrixSigma(freqs,opts)
%% tmatrixSigma - medium cross-sections from the sphere T-matrix
%
%   out = tmatrixSigma([300e3 500e3 800e3])
%
% Builds Sigma_ij(theta) = n0 * |f_ij(theta)|^2 * S(q_ij(theta)) and from it
% the transport mean free paths and the diffusivity, with the same formulas
% MaterialClass.prepareSigma uses for the Born path, so the two are compared
% like with like.
%
% q depends on the mode pair:
%   q_PP = 2 k_P sin(th/2)      q_SS = 2 k_S sin(th/2)
%   q_PS = sqrt(k_P^2 + k_S^2 - 2 k_P k_S cos th)
%
% S(q) comes from the sphere centres (packing), not from theory.

if nargin < 1, freqs = [300e3 500e3 800e3]; end
if nargin < 2, opts = struct; end
if ~isfield(opts,'phi'),   opts.phi   = 0.4003484; end   % from the covariance file
if ~isfield(opts,'a'),     opts.a     = 4e-3;      end
if ~isfield(opts,'ndir'),  opts.ndir  = 200;       end
if ~isfield(opts,'nq'),    opts.nq    = 80;        end
if ~isfield(opts,'centres')
    % the ten realizations behind the autocorrelation used by the RTE
    % (the RTE code uses mean(autoco,2), the average over the ten)
    raw = '/Users/lac/Nextcloud/ArtigoManda/DescricaoEstatistica/Raw/';
    opts.centres = arrayfun(@(i) sprintf('%ssphere_centers_45P_LSMETHOD_R%d.dat',raw,i), ...
                            1:10,'UniformOutput',false);
end

matriz = struct('lam',1.1229e10,'mu',1.0378e10,'rho',2050);
nucleo = struct('lam',2.5774e10,'mu',2.8322e10,'rho',2450);

a  = opts.a;
n0 = opts.phi/((4/3)*pi*a^3);
vp = sqrt((matriz.lam+2*matriz.mu)/matriz.rho);
vs = sqrt(matriz.mu/matriz.rho);

fprintf('n0 = %.4g /m3 (phi = %.4f, a = %.1f mm) | vp=%.0f vs=%.0f\n', ...
        n0,opts.phi,a*1e3,vp,vs);

%% structure factor of the centres
qmax = 2*2*pi*max(freqs)/vs*1.05;
qmin = 2*pi/ (0.15);                 % larger than 2pi/L of the box
[qg,Sg] = TMatrix.structureFactor(opts.centres,qmin,qmax,opts.nq,opts.ndir);
Sfun = @(q) interp1(qg,Sg,min(max(q,qg(1)),qg(end)),'linear','extrap');

fprintf('S(q): %d points from %.0f to %.0f 1/m | S at the peak = %.2f\n\n', ...
        numel(qg),qg(1),qg(end),max(Sg));

%% per band
out = struct('f',{},'Sigma',{},'Sigmapr',{},'lstar',{},'D',{},'P2P',{},'S2S',{});
fprintf('%8s %10s %10s %10s %10s %10s\n', ...
        'f (kHz)','l*_P (mm)','l*_S (mm)','S*_S (1/m)','D (m2/s)','P2P');
for i = 1:numel(freqs)
    om = 2*pi*freqs(i);
    sc = TMatrix.sphereScattering(a,matriz,nucleo,om);

    kp = om/vp; ks = om/vs;
    th = sc.th;

    qPP = 2*kp*sin(th/2);
    qSS = 2*ks*sin(th/2);
    qPS = sqrt(kp^2 + ks^2 - 2*kp*ks*cos(th));

    sg = cell(2);
    sg{1,1} = n0*sc.dsig_PP.*Sfun(qPP);
    sg{1,2} = n0*sc.dsig_PS.*Sfun(qPS);
    sg{2,1} = n0*sc.dsig_SP.*Sfun(qPS);
    sg{2,2} = n0*sc.dsig_SS.*Sfun(qSS);

    Sig = zeros(2); Spr = zeros(2);
    for r = 1:2
        for c = 1:2
            Sig(r,c) = trapz(th,2*pi*sin(th).*sg{r,c});
            Spr(r,c) = trapz(th,2*pi*sin(th).*cos(th).*sg{r,c});
        end
    end

    % prepareSigma works with Sigma in 1/s (collision rate), not in 1/m: there
    % Sigma comes from the Born formula, which already carries omega. Here
    % Sigma is in 1/m, so it is multiplied by the velocity of the INCIDENT
    % mode (the row).
    Sig = [vp;vs].*Sig;
    Spr = [vp;vs].*Spr;

    % same formulas as prepareSigma
    den = (Sig(1,1)+Sig(1,2)-Spr(1,1))*(Sig(2,2)+Sig(2,1)-Spr(2,2)) - Spr(1,2)*Spr(2,1);
    lP  = (vp*(Sig(2,2)+Sig(2,1)-Spr(2,2)) + vs*Spr(1,2))/den;
    lS  = (vs*(Sig(1,1)+Sig(1,2)-Spr(1,1)) + vp*Spr(2,1))/den;
    K   = vp/vs;
    D   = (vp*lP/3 + 2*K^3*vs*lS/3)/(1+2*K^3);

    o = struct('f',freqs(i),'Sigma',Sig,'Sigmapr',Spr, ...
               'lstar',[lP lS],'D',D, ...
               'P2P',Sig(1,1)/sum(Sig(1,:)),'S2S',Sig(2,2)/sum(Sig(2,:)));
    out(end+1) = o; %#ok<AGROW>

    fprintf('%8.0f %10.2f %10.2f %10.2f %10.2f %10.3f\n', ...
            freqs(i)/1e3,lP*1e3,lS*1e3,1/lS,D,o.P2P);
end

fprintf(['\nmeasured reference:  Sigma*_SEM = 15.1 / 15.8 / 17.6 1/m  (+17%% from LF to HF)\n' ...
         '                    Sigma*_Born = 18.2 / 22.6 / 30.8 1/m  (+69%%)\n' ...
         '                    D_SEM = 65.6 / 62.4 / 56.2 | D_Born = 54.2 / 43.7 / 32.1 m2/s\n']);
end
