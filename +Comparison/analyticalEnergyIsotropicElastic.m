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

    kk = [0 besselroots(0,100).'/r]; Nk = length(kk);
    Nw = 2^12; Lw = 100; ww = linspace(-Lw/2,Lw/2,Nw); dww = Lw/Nw;
    dt = 2*pi/Lw; tt = (-Nw/2:Nw/2-1)*dt;

    out1 = zeros(Nw,Nk-1); 
    fun = @(x,w) x.*besselj(0,x.*r).*Emp(x,-1i*w,K);
    parfor j=1:Nk-1
        out1(:,j) = integral(@(x) fun(x,ww'),kk(j),kk(j+1),'RelTol',0,'AbsTol',1e-10,'ArrayValued',true);
    end
    temp = forwardAverage(cumsum(out1,2),Nk-2);
    
    Epm = 1/(2*pi)^2*interp1(tt,abs(dww*fftshift(fft(temp(:,end)))),a);
    Epm(~H(a-r)) = 0;
    Ep = Ep0 + Ep1 + Epm;

    out1 = zeros(Nw,Nk-1); 
    fun = @(x,w) x.*besselj(0,x.*r).*Ems(x,-1i*w,K);
    parfor j=1:Nk-1
        out1(:,j) = integral(@(x) fun(x,ww'),kk(j),kk(j+1),'RelTol',0,'AbsTol',1e-10,'ArrayValued',true);
    end
    temp = forwardAverage(cumsum(out1,2),Nk-2);

    Esm = 1/(2*pi)^2*interp1(tt,abs(dww*fftshift(fft(temp(:,end)))),a);
    Esm(~H(a-r)) = 0;
    Es = Es0 + Es1 + Esm;

    E = zeros(length(time),2);
    E(:,1) = Ep; E(:,2) = Es;

elseif dim == 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Direct waves (no scattering)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gaussian pulse used to define a Dirac delta function
    E0 = (deltaFunction(a-r',width) + 1.5*K^5*deltaFunction(a-K*r',width)) ...
                               .*exp(-a).*1./(4*pi*r'.^2*(1+1.5*K^5));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Single scattering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Epp1 = 1*Sigmapp*K_func(a./r');
    Epp1(~H(a-r')) = 0;

    Eps1 = 1*K*Sigmaps*Kc_func(a./r',K);
    Eps1(a./r'==1) = 1*K*Sigmaps*atanh(1/K);
    Eps1(~H(a-r')) = 0;

    Esp1 = 1.5*K^5*K*Sigmasp*Kc_func(a./r',K);
    Esp1(a./r'==1) = 1.5*K^5*K*Sigmasp*atanh(1/K);
    Esp1(~H(a-r')) = 0;

    Ess1 = 1.5*K^5*K*Sigmass*K_func(a/K./r');
    Ess1(~H(a-K*r')) = 0;

    E1 = (Epp1+Eps1+Esp1+Ess1).*exp(-a)./(4*pi*r'.^2*eta*(1+1.5*K^5));

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
    Nk = 2^14; Nw = 2^14;
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
    [Xq,Yq] = meshgrid(r,a);
    E_interp = interp2(rr,tt,E_fourier,Xq,Yq);

    E_multiple = eta./((1+1.5*K^5)*(2*pi)^2*r').*E_interp';
    E_multiple(~H(a-r')) = 0;

    E = E0 + E1 + E_multiple;

else
    error('Dimension should be either 2 or 3')
end
end

function out = M_func(x,K)
% Computes the M function introduced in Nakahara & Yoshimito 2011
out = zeros(size(x));
out(x<1) = 0;
ind1 = x>=K;
out(ind1) = ellipke(sqrt((K^2-1)./(x(ind1).^2-1)))./sqrt(x(ind1).^2-1);
ind2 = x>=1 & x<K;
out(ind2) = ellipke(sqrt((x(ind2).^2-1)/(K.^2-1)))/sqrt(K^2-1);
end

function output = forwardAverage(input, m)
% Computes the forward repeated averaging operator (Euler transformation)
    currentVector = input;    
    % Apply the forward averaging m times
    for i = 1:m
        % Check if the current vector has enough elements
        if size(currentVector,2) < 2
            break;
        end
        % Compute the forward average
        currentVector = (currentVector(:,1:end-1) + currentVector(:,2:end)) / 2;
    end    
    output = currentVector;
end

function j = besselroots(v, n) 
% Trivial case:
if ( n == 0 )
    j = [];
    return
end
% Check inputs:
if ( ~isscalar(n) || (round(n) - n ~= 0) || n < 0 )
    error('CHEBFUN:besselroots:inputN', 'Input N must be a positive integer');
end
% McMahon's expansion. This expansion gives very accurate approximation 
% for the sth zero (s >= 7) in the whole region V >=- 1, and moderate
% approximation in other cases.
s = (1:n)';
mu = 4*v^2;
a1 = 1 / 8;
a3 = (7*mu-31) / 384;
a5 = 4*(3779+mu*(-982+83*mu)) / 61440; % Evaluate via Horner's method.
a7 = 6*(-6277237+mu*(1585743+mu*(-153855+6949*mu))) / 20643840;
a9 = 144*(2092163573+mu*(-512062548+mu*(48010494+mu*(-2479316+70197*mu)))) ...
     / 11890851840;
a11 = 720*(-8249725736393+mu*(1982611456181+mu*(-179289628602+mu*(8903961290 + ...
    mu*(-287149133+5592657*mu))))) / 10463949619200;
a13 = 576*(423748443625564327 + mu*(-100847472093088506+mu*(8929489333108377 + ...
    mu*(-426353946885548+mu*(13172003634537+mu*(-291245357370 + mu*4148944183)))))) ...
     / 13059009124761600;
b = .25*(2*v+4*s-1)*pi; % beta
j = b - (mu-1)*polyval([a13 0 a11 0 a9 0 a7 0 a5 0 a3 0 a1 0], 1./b);
if ( v == 0 )
    % First 20 roots of J0(x) are precomputed (using Wolfram Alpha):
    j(1:20) = [
        2.4048255576957728
        5.5200781102863106
        8.6537279129110122
        11.791534439014281
        14.930917708487785
        18.071063967910922
        21.211636629879258
        24.352471530749302
        27.493479132040254
        30.634606468431975
        33.775820213573568
        36.917098353664044
        40.058425764628239
        43.199791713176730
        46.341188371661814
        49.482609897397817
        52.624051841114996
        55.765510755019979
        58.906983926080942
        62.048469190227170];
elseif ( v >= -1 && v <= 5 )
    % Piessens's Chebyshev series approximations (1984). Calculates the 6 first
    % zeros to at least 12 decimal figures in region -1 <= V <= 5:
    C = [
       2.883975316228  8.263194332307 11.493871452173 14.689036505931 17.866882871378 21.034784308088
       0.767665211539  4.209200330779  4.317988625384  4.387437455306  4.435717974422  4.471319438161
      -0.086538804759 -0.164644722483 -0.130667664397 -0.109469595763 -0.094492317231 -0.083234240394
       0.020433979038  0.039764618826  0.023009510531  0.015359574754  0.011070071951  0.008388073020
      -0.006103761347 -0.011799527177 -0.004987164201 -0.002655024938 -0.001598668225 -0.001042443435
       0.002046841322  0.003893555229  0.001204453026  0.000511852711  0.000257620149  0.000144611721
      -0.000734476579 -0.001369989689 -0.000310786051 -0.000105522473 -0.000044416219 -0.000021469973
       0.000275336751  0.000503054700  0.000083834770  0.000022761626  0.000008016197  0.000003337753
      -0.000106375704 -0.000190381770 -0.000023343325 -0.000005071979 -0.000001495224 -0.000000536428
       0.000042003336  0.000073681222  0.000006655551  0.000001158094  0.000000285903  0.000000088402
      -0.000016858623 -0.000029010830 -0.000001932603 -0.000000269480 -0.000000055734 -0.000000014856
       0.000006852440  0.000011579131  0.000000569367  0.000000063657  0.000000011033  0.000000002536
      -0.000002813300 -0.000004672877 -0.000000169722 -0.000000015222 -0.000000002212 -0.000000000438
       0.000001164419  0.000001903082  0.000000051084  0.000000003677  0.000000000448  0.000000000077
      -0.000000485189 -0.000000781030 -0.000000015501 -0.000000000896 -0.000000000092 -0.000000000014
       0.000000203309  0.000000322648  0.000000004736  0.000000000220  0.000000000019  0.000000000002
      -0.000000085602 -0.000000134047 -0.000000001456 -0.000000000054 -0.000000000004               0
       0.000000036192  0.000000055969  0.000000000450  0.000000000013               0               0
      -0.000000015357 -0.000000023472 -0.000000000140 -0.000000000003               0               0
       0.000000006537  0.000000009882  0.000000000043  0.000000000001               0               0
      -0.000000002791 -0.000000004175 -0.000000000014               0               0               0
       0.000000001194  0.000000001770  0.000000000004               0               0               0
      -0.000000000512 -0.000000000752               0               0               0               0
       0.000000000220  0.000000000321               0               0               0               0
      -0.000000000095 -0.000000000137               0               0               0               0
       0.000000000041  0.000000000059               0               0               0               0
      -0.000000000018 -0.000000000025               0               0               0               0
       0.000000000008  0.000000000011               0               0               0               0
      -0.000000000003 -0.000000000005               0               0               0               0
       0.000000000001  0.000000000002               0               0               0               0];
    j(1:6) = chebtech.clenshaw((v-2)/3, C).'; % Evaluate via Clenshaw.
    j(1) = j(1) * sqrt(v+1);                  % Scale the first root.
end
j = j(1:n);                                   % Trim unnecessary points.
end