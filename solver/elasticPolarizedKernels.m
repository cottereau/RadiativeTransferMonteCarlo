function kernel = elasticPolarizedKernels(mat,n,m,incident,outgoing)
%ELASTICPOLARIZEDKERNELS 3D isotropic Born P/S scattering Kraus map.
% n,m are incident/outgoing unit directions. incident/outgoing are 'P'/'S'.
% For an incident coherency W, outgoing coherency per solid angle is
% sum_l K_l W K_l'. Rates have units 1/(s sr). S bases are deterministic
% transverse frames; changing either frame conjugates the map covariantly.
% The real covariance ensemble is represented by its eigenfields, so a
% pure-state Monte Carlo unraveling preserves the full mixed-state map.
% Coefficients and PSD convention match MaterialClass.calcSigmaIsotropic.
if ~isa(mat,'MaterialClass') || ~isequal(mat.acoustics,false) || mat.d~=3 || ...
        ~strcmp(mat.scatteringModel,'parameterFluctuations') || ~strcmp(mat.correlationStructure,'isotropic') || ~isa(mat.Phi,'function_handle')
    error('Fonte:PolarizedMaterial','A 3D isotropic elastic Born material with Phi is required.');
end
validateattributes(n,{'numeric'},{'real','finite','numel',3});
validateattributes(m,{'numeric'},{'real','finite','numel',3});
n=double(n(:)); m=double(m(:));
if abs(norm(n)-1)>1e-10 || abs(norm(m)-1)>1e-10
    error('Fonte:PolarizedDirections','Directions must be unit vectors.');
end
if ~any(strcmp(incident,{'P','S'})) || ~any(strcmp(outgoing,{'P','S'}))
    error('Fonte:PolarizedModes','Modes must be P or S.');
end
validateattributes(mat.vp,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(mat.vs,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(mat.Frequency,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(mat.CorrelationLength,{'numeric'},{'real','finite','scalar','positive'});
cv=double(mat.coefficients_of_variation(:)); corr=double(mat.correlation_coefficients(:));
if numel(cv)~=3 || numel(corr)~=3 || any(~isfinite([cv;corr])) || any(cv<0)
    error('Fonte:PolarizedCovariance','Need three finite nonnegative CVs and three correlations.');
end
R=[1 corr(1) corr(2);corr(1) 1 corr(3);corr(2) corr(3) 1];
C=(cv*cv').*R; C=(C+C')/2;
[U,D]=eig(C); eigval=diag(D);
if min(eigval)<-1e-12*max(1,max(eigval))
    error('Fonte:PolarizedCovariance','Fluctuation covariance is not positive semidefinite.');
end
mix=U*diag(sqrt(max(eigval,0)));
Kratio=double(mat.vp)/double(mat.vs); c=n'*m;
omega=2*pi*double(mat.Frequency); zp=omega*double(mat.CorrelationLength)/double(mat.vp);
zs=omega*double(mat.CorrelationLength)/double(mat.vs);
if strcmp(incident,'P'), zi=zp; else, zi=zs; end
if strcmp(outgoing,'P'), zo=zp; else, zo=zs; end
q=sqrt(max(0,zi^2+zo^2-2*zi*zo*c));
phi=mat.Phi(q);
if ~isscalar(phi) || ~isreal(phi) || ~isfinite(phi) || phi<0
    error('Fonte:PolarizedPSD','Phi must be a finite nonnegative real spectral density.');
end
if strcmp(incident,'P'), Bin=n; else, Bin=transverseBasis(n); end
if strcmp(outgoing,'P'), Bout=m; else, Bout=transverseBasis(m); end
di=size(Bin,2); do=size(Bout,2);
A=zeros(do,di,3); base=(pi/2)*omega*zp^3;
if strcmp(incident,'P') && strcmp(outgoing,'P')
    A(1,1,1)=1-2/Kratio^2;
    A(1,1,2)=2*c^2/Kratio^2;
    A(1,1,3)=-c;
    pref=base;
elseif strcmp(incident,'P') && strcmp(outgoing,'S')
    projection=Bout'*n;
    A(:,:,2)=2*c*projection;
    A(:,:,3)=-Kratio*projection;
    pref=base*Kratio;
elseif strcmp(incident,'S') && strcmp(outgoing,'P')
    projection=m'*Bin;
    A(:,:,2)=(2*c/Kratio)*projection;
    A(:,:,3)=-projection;
    pref=base;
else
    for a=1:do
        f=Bout(:,a);
        for b=1:di
            e=Bin(:,b);
            A(a,b,2)=(f'*n)*(m'*e)+c*(f'*e);
            A(a,b,3)=-(f'*e);
        end
    end
    pref=(pi/2)*omega*zs^3;
end
maps=zeros(do,di,3);
for ell=1:3
    for a=1:3
        maps(:,:,ell)=maps(:,:,ell)+mix(a,ell)*A(:,:,a);
    end
end
maps=sqrt(pref*phi)*maps;
rate=zeros(di);
for ell=1:3, rate=rate+maps(:,:,ell)'*maps(:,:,ell); end
kernel=struct('kraus',maps,'rateMatrix',rate, ...
    'unpolarizedRate',trace(rate)/di,'inputBasis',Bin, ...
    'outputBasis',Bout,'incident',incident,'outgoing',outgoing, ...
    'qCorrelation',q,'units','1/(s sr)');
end

function B=transverseBasis(n)
[~,j]=min(abs(n)); ref=zeros(3,1); ref(j)=1;
u=cross(n,ref); u=u/norm(u); B=[u cross(n,u)];
end
