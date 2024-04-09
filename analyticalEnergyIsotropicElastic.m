function E = analyticalEnergyIsotropicElastic(dim,K,r,time,Sigma)
% This function computes the normalized P and S-wave energy densities for
% 2D and 3D elastic isotropic scattering media, based on H. Nakahara and K.
% Yoshimoto, 2011 (for 2D case) and H. Sato 1994 (for 3D case)

% Inputs
% dim   : dimension of the problem (2 or 3)
% K     : vp/vs (ratio between P and S-wave velocities)
% r     : normalized source-station distance
% Sigma : cell array containig total scattering cross sections

% Output
% E : normalized (by (eta/vp)^2) energy density
%     if d=2, an array where the second dimension specifies the energy of
%     each polarization type (first dimension = 'P', second dimension = 'S')

% Assumption : Following Sato 1994, this code assumes Sigmap = Sigmas = eta

Sigmapp = Sigma{1,1}; Sigmaps = Sigma{1,2};
Sigmasp = Sigma{2,1}; Sigmass = Sigma{2,2};
Sigmap = Sigmapp + Sigmaps;
eta = Sigmap;

% normalized time (a = eta*t)
a = eta*time;

% Some functions which are subsequently used
% Gaussian pulse used to define a Dirac delta function
width = 0.01; % Width of the Gaussian pulse
deltaFunction = @(x,width) exp(-x.^2/(2*width^2))/(sqrt(2*pi)*width);
H = @(x) x>0; % Modified Heviside function
Kc_func = @(x,K) atanh(x/K)./x.*(H(x-1)-H(x-K)) + ...
    (atanh(1./x)+atanh(K./x))./x.*H(x-K);
K_func = @(x) 2*atanh(1./x)./x.*H(x-1);

if dim==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Direct waves (no scattering)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ep0 = exp(-a)./(2*pi*r).*deltaFunction(a-r,width);
    Es0 = 1.5*K^5*K*exp(-a)./(2*pi*r).*deltaFunction(a-K*r,width);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Single scattering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Epp1 = (Sigmapp/eta)*exp(-a)./(2*pi*r.*sqrt((a./r).^2 - 1));
    Epp1(~H(a-r)) = 0;
    Esp1 = 1.5*K^5*K*(Sigmasp/eta)*exp(-a).*M_func(a./r,K)./(pi^2*r);
    Esp1(~H(a-r)) = 0;

    Ep1 = Epp1 + Esp1;

    Ess1 = (Sigmass/eta)*1.5*K^5*K*exp(-a)./(2*pi*r.*sqrt((a./K/r).^2 - 1));
    Ess1(~H(a-K*r)) = 0;
    Eps1 = 1*K*(Sigmaps/eta)*exp(-a).*M_func(a./r,K)./(pi^2*r);
    Eps1(~H(a-r)) = 0;

    Es1 = Ess1 + Eps1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multiple scattering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Gp = @(k,s) (1/eta)*1./sqrt((1+s).^2+k.^2);
    Gs = @(k,s,K) (1/eta)*1./sqrt((1+s).^2+(k/K).^2);

    Emp = @(k,s,K) ( 1*Gp(k,s).*( (Sigmapp*Gp(k,s)).^2.*(1-Sigmass*Gs(k,s,K)) + ...
        Sigmaps*Sigmasp*(1+Sigmapp*Gp(k,s)).*Gp(k,s).*Gs(k,s,K) ) ...
        + 1.5*K^5*Gs(k,s,K).*( Sigmasp*Gp(k,s).*( Sigmapp*Gp(k,s).*...
        (1-Sigmass*Gs(k,s,K))+Sigmass*Gs(k,s,K)+Sigmaps*Sigmasp*Gp(k,s).*Gs(k,s,K) ) ) ) ...
        ./ ( (1-Sigmapp*Gp(k,s)).*(1-Sigmass*Gs(k,s,K))-Sigmaps*Sigmasp*Gp(k,s).*Gs(k,s,K) );

    Ems = @(k,s,K) ( 1.5*K^5*Gs(k,s,K).*( (Sigmass*Gs(k,s,K)).^2.*(1-Sigmapp*Gp(k,s)) + ...
        Sigmaps*Sigmasp*(1+Sigmass*Gs(k,s,K)).*Gp(k,s).*Gs(k,s,K) ) ...
        + 1*Gp(k,s).*( Sigmaps*Gs(k,s,K).*( Sigmass*Gs(k,s,K).*...
        (1-Sigmapp*Gp(k,s))+Sigmapp*Gp(k,s)+Sigmaps*Sigmasp*Gp(k,s).*Gs(k,s,K) ) ) ) ...
        ./ ( (1-Sigmapp*Gp(k,s)).*(1-Sigmass*Gs(k,s,K))-Sigmaps*Sigmasp*Gp(k,s).*Gs(k,s,K) );

    Lk = 500; Lw = 500;
    Nk = 2^12; Nw = 2^12; dww = Lw/Nw;
    kk = linspace(0,Lk,Nk);
    ww = linspace(-Lw/2,Lw/2,Nw);

    dt = 2*pi/Lw;
    tt = (-Nw/2:Nw/2-1)*dt;

    E_fourier_p = abs(dww*fftshift(fft(Emp(kk',-1i*ww,K),[],2),2));
    E_fourier_s = abs(dww*fftshift(fft(Ems(kk',-1i*ww,K),[],2),2));
    Epm_interp = interp1(tt,E_fourier_p',a);
    Esm_interp = interp1(tt,E_fourier_s',a);
    Epm = eta/(2*pi)^2*trapz(kk,kk.*besselj(0,kk.*r).*Epm_interp,2);
    Epm(~H(a-r)) = 0;

    Esm = eta/(2*pi)^2*trapz(kk,kk.*besselj(0,kk.*r).*Esm_interp,2);
    Esm(~H(a-r)) = 0;

    Ep = Ep0 + Ep1 + Epm';
    Es = Es0 + Es1 + Esm';
    E = zeros(length(time),2);
    E(:,1) = Ep; E(:,2) = Es;

elseif dim == 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Direct waves (no scattering)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gaussian pulse used to define a Dirac delta function
    width = 0.01; % Width of the Gaussian pulse
    deltaFunction = @(x) exp(-x.^2/(2*width^2))/(sqrt(2*pi)*width);
    E0 = (deltaFunction(a-r) + 1.5*K^5*deltaFunction(a-K*r)) ...
                               .*exp(-a).*1./(4*pi*r^2*(1+1.5*K^5));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Single scattering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Epp1 = 1*Sigmapp*K_func(a/r);
    Epp1(~H(a-r)) = 0;

    Eps1 = 1*K*Sigmaps*Kc_func(a/r,K);
    Eps1(a/r==1) = 1*K*Sigmaps*atanh(1/K);
    Eps1(~H(a-r)) = 0;

    Esp1 = 1.5*K^5*K*Sigmasp*Kc_func(a/r,K);
    Esp1(a/r==1) = 1.5*K^5*K*Sigmasp*atanh(1/K);
    Esp1(~H(a-r)) = 0;

    Ess1 = 1.5*K^5*K*Sigmass*K_func(a/K/r);
    Ess1(~H(a-K*r)) = 0;

    E1 = (Epp1+Eps1+Esp1+Ess1).*exp(-a)./(4*pi*r^2*eta*(1+1.5*K^5));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Multiple scattering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Gp = @(k,s) (1/eta./k).*atan(k./(s+1));
    Gs = @(k,s,K) (K/eta./k).*atan(k./(K*(s+1)));

    Em = @(k,s,K) (1*Gp(k,s).*(Sigmaps*Gs(k,s,K).*(Sigmass*Gs(k,s,K).*(1-Sigmapp*Gp(k,s))+Sigmasp*Gp(k,s).*(1+Sigmaps*Gs(k,s,K))) + ...
        Sigmapp*Gp(k,s).*(Sigmapp*Gp(k,s).*(1-Sigmass*Gs(k,s,K))+Sigmaps*Gs(k,s,K).*(1+Sigmasp*Gp(k,s)))) + ...
        +1.5*K^5*Gs(k,s,K).*(Sigmasp*Gp(k,s).*(Sigmapp*Gp(k,s).*(1-Sigmass*Gs(k,s,K))+Sigmaps*Gs(k,s,K).*(1+Sigmasp*Gp(k,s))) + ...
        Sigmass*Gs(k,s,K).*(Sigmass*Gs(k,s,K).*(1-Sigmapp*Gp(k,s))+Sigmasp*Gp(k,s).*(1+Sigmaps*Gs(k,s,K)))))./ ...
        ( (1-Sigmapp*Gp(k,s)).*(1-Sigmass*Gs(k,s,K))-Sigmaps*Sigmasp*Gp(k,s).*Gs(k,s,K) );
    Em0 = @(s,K) (1*(1/eta./(s+1)).*(Sigmaps*(1/eta./(s+1)).*(Sigmass*(1/eta./(s+1)).*(1-Sigmapp*(1/eta./(s+1)))+Sigmasp*(1/eta./(s+1)).*(1+Sigmaps*(1/eta./(s+1)))) + ...
        Sigmapp*(1/eta./(s+1)).*(Sigmapp*(1/eta./(s+1)).*(1-Sigmass*(1/eta./(s+1)))+Sigmaps*(1/eta./(s+1)).*(1+Sigmasp*(1/eta./(s+1))))) + ...
        +1.5*K^5*(1/eta./(s+1)).*(Sigmasp*(1/eta./(s+1)).*(Sigmapp*(1/eta./(s+1)).*(1-Sigmass*(1/eta./(s+1)))+Sigmaps*(1/eta./(s+1)).*(1+Sigmasp*(1/eta./(s+1)))) + ...
        Sigmass*(1/eta./(s+1)).*(Sigmass*(1/eta./(s+1)).*(1-Sigmapp*(1/eta./(s+1)))+Sigmasp*(1/eta./(s+1)).*(1+Sigmaps*(1/eta./(s+1))))))./ ...
        ( (1-Sigmapp*(1/eta./(s+1))).*(1-Sigmass*(1/eta./(s+1)))-Sigmaps*Sigmasp*(1/eta./(s+1)).*(1/eta./(s+1)) );

    % Definition of normalized variables k and w
    Lk = 500; Lw = 500;
    Nk = 2^12; Nw = 2^12;
    dkk = Lk/Nk; dww = Lw/Nw;
    kk = linspace(-Lk/2,Lk/2,Nk);
    ww = linspace(-Lw/2,Lw/2,Nw);

    y = Em(-kk',-1i*ww,K);
    if any(kk'==0)
        y(kk'==0,:) = Em0(-1i*ww,K);
    end

    dr = 2*pi/Lk; dt = 2*pi/Lw;
    rr = (-Nk/2:Nk/2-1)*dr;
    tt = (-Nw/2:Nw/2-1)*dt;

    E_fourier = abs( dkk*dww*fftshift( fft2( 1i*kk'/(2*pi).*y ) ) );
    E_interp = interp2(rr,tt,E_fourier,r,a);

    E_multiple = eta/((1+1.5*K^5)*(2*pi)^2*r).*E_interp';
    E_multiple(~H(a-r)) = 0;

    E = E0 + E1 + E_multiple;

else
    error('Dimension should be either 2 or 3')
end


function out = M_func(x,K)

out = zeros(size(x));
out(x<1) = 0;
ind1 = x>=K;
out(ind1) = ellipke(sqrt((K^2-1)./(x(ind1).^2-1)))./sqrt(x(ind1).^2-1);
ind2 = x>=1 & x<K;
out(ind2) = ellipke(sqrt((x(ind2).^2-1)/(K.^2-1)))/sqrt(K^2-1);
