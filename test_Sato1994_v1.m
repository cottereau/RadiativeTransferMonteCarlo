% Computation of normalized energy densities for 3D elastic isotropic 
% scattering media, based on Sato 1994

material = struct( 'acoustics', false, ...
                   'vp', 6, ...
                   'vs', 3.46);
vp = material.vp; vs = material.vs;
gama = vp/vs;

% Numerical case study in Sato 1994 (Fig. 4, left)
eta = 0.1;
Sigmapp = eta/2; Sigmaps = eta/2;
Sigmasp = eta/2; Sigmass = eta/2;

H = @(x) x>0; % Modified Heviside function
Kc_func = @(x,K) atanh(x/K)./x.*(H(x-1)-H(x-K)) + (atanh(1./x)+atanh(K./x))./x.*H(x-K);
K_func = @(x) 2*atanh(1./x)./x.*H(x-1);

% normalized time (a = eta*t)
Nt = 256;
a = linspace(0,4,Nt);
[m,n] = size(a); if m>n, a=a'; end

% normalized distance (b = sigma*r/vp)
Nr = 256;
b = linspace(0,2,Nr);
[m,n] = size(b); if m<n, b=b'; end

% Direct waves (no scattering)
E0 = zeros(Nr,Nt);
[~,ind1] = min((a-b).^2,[],1);
for i = 1:Nt
    E0(ind1(i),i) = E0(ind1(i),i) + exp(-a(i)).*1./(4*pi*b(i).^2*(1+1.5*gama^2));
end

[~,ind2] = min((a-gama*b).^2,[],1);
for i = 1:Nt
    E0(ind2(i),i) = E0(ind2(i),i) + exp(-a(i)).*1.5*gama^5./(4*pi*b(i).^2*(1+1.5*gama^2));
end

% Single scattering terms
E1_ps = 1*gama*Sigmaps*Kc_func(a./b,gama);
E1_ps(H(a-b)) = 0;

E1_pp = 1*Sigmapp*K_func(a./b);
E1_pp(H(a-b)) = 0;

E1_sp = 1.5*gama^5*gama*Sigmasp*Kc_func(a./b,gama);
E1_sp(H(a-b)) = 0;

E1_ss = 1.5*gama^5*gama*Sigmass*K_func(a/gama./b);
E1_ss(H(a-gama*b)) = 0;

E1 = (E1_ps+E1_pp+E1_sp+E1_ss).*exp(-a)./(4*pi*b.^2*eta*(1+1.5*gama^5));

% Multiple scattering term
Nk = 2^9;
kk = linspace(-5,5,Nk); 
dkk = mean(diff(kk)); fsk = 1/dkk;
rr = (-Nk/2:Nk/2-1)*fsk/Nk;

Nw = 2^8;
ww = linspace(-5,5,Nw); 
dww = mean(diff(ww)); fsw = 1/dww; 
tt = (-Nw/2:Nw/2-1)*fsw/Nw;

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
y = Em(-kk.',-1i*ww,gama);
if any(kk==0)
    y(kk==0,:) = Em0(-1i*ww,gama);
end

Y = dkk*dww*fft2(1i*kk.'.*y/2/pi);

E_multiple = (eta/(1+1.5*gama^5)/(2*pi)^2./rr.').*abs(Y);

[X,Y] = meshgrid(tt,rr);
Em_final = interp2(X,Y,E_multiple,a,b);

% superposition of contributions from different scattering orders
E = E0 + E1 + Em_final;

r_target = 1;
inds = find(abs(b-r_target)<0.01);
figure; plot(a,abs(E(inds(1),:)))