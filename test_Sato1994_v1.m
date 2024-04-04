% Computation of normalized energy densities for 3D elastic isotropic 
% scattering media, based on Sato 1994

gama = sqrt(3);

% Numerical case study in Sato 1994 (Fig. 4, left)
eta = 1;
Sigmapp = eta/2; Sigmaps = eta/2;
Sigmasp = eta/2; Sigmass = eta/2;

H = @(x) x>0; % Modified Heviside function
Kc_func = @(x,K) atanh(x/K)./x.*(H(x-1)-H(x-K)) + (atanh(1./x)+atanh(K./x))./x.*H(x-K);
K_func = @(x) 2*atanh(1./x)./x.*H(x-1);

% normalized time (a = eta*t)
Nt = 2^10;
a = linspace(0.01,8,Nt);
[m,n] = size(a); if m>n, a=a'; end

% normalized distance (b = sigma*r/vp)
Nr = 2^10;
b = linspace(0.01,8,Nr);
[m,n] = size(b); if m<n, b=b'; end

[B,A] = meshgrid(b,a);

% Direct waves (no scattering)
% Gaussian pulse used to define a Dirac delta function
width = 0.01; % Width of the Gaussian pulse
deltaFunction_1 = exp(-(a-b).^2/(2*width^2))/(sqrt(2*pi)*width);
deltaFunction_2 = exp(-(a-gama*b).^2/(2*width^2))/(sqrt(2*pi)*width);
E0 = (deltaFunction_1 + 1.5*gama^5*deltaFunction_2).*exp(-a).*1./(4*pi*b.^2*(1+1.5*gama^5));

% E0 = zeros(Nr,Nt);
% [~,ind1] = min((a-b).^2,[],1);
% for i = 1:Nt
%     E0(ind1(i),i) = E0(ind1(i),i) + exp(-a(i)).*1./(4*pi*b(i).^2*(1+1.5*gama^5));
% end

% [~,ind2] = min((a-gama*b).^2,[],1);
% for i = 1:Nt
%     E0(ind2(i),i) = E0(ind2(i),i) + exp(-a(i)).*1.5*gama^5./(4*pi*b(i).^2*(1+1.5*gama^5));
% end

% Single scattering terms
E1_pp = 1*Sigmapp*K_func(a./b);
E1_pp(~H(a-b)) = 0;

E1_ps = 1*gama*Sigmaps*Kc_func(a./b,gama);
E1_ps(a./b==1) = 1*gama*Sigmaps*atanh(1/gama);
E1_ps(~H(a-b)) = 0;

E1_sp = 1.5*gama^5*gama*Sigmasp*Kc_func(a./b,gama);
E1_sp(a./b==1) = 1.5*gama^5*gama*Sigmasp*atanh(1/gama);
E1_sp(~H(a-b)) = 0;

E1_ss = 1.5*gama^5*gama*Sigmass*K_func(a/gama./b);
E1_ss(~H(a-gama*b)) = 0;

E1 = (E1_pp+E1_ps+E1_sp+E1_ss).*exp(-a)./(4*pi*b.^2*eta*(1+1.5*gama^5));

% Multiple scattering term
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

y = Em(-kk',-1i*ww,gama);
if any(kk'==0)
    y(kk'==0,:) = Em0(-1i*ww,gama);
end

dr = 2*pi/Lk; dt = 2*pi/Lw;
rr = (-Nk/2:Nk/2-1)*dr;
tt = (-Nw/2:Nw/2-1)*dt;
[R,T] = meshgrid(rr,tt);

E_fourier = abs( dkk*dww*fftshift( fft2( 1i*kk'/(2*pi).*y ) ) );
E_interp = interp2(R,T,E_fourier,B,A);
E_multiple = eta./((1+1.5*gama^5)*(2*pi)^2*b).*E_interp;
E = E0 + E1 + E_multiple;

r_target = 0.5;
inds = find(abs(b-r_target)<0.005);
figure; plot(a,E(inds(1),:)); xlim([0 4]);
