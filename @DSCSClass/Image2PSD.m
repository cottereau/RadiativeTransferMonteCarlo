function [out,normk] = Image2PSD(Im,dx,dy,theta)

% check dx and dy
if ~exist('dx','var')
    dx = 1;
end
if ~exist('dy','var')
    dy = 1;
end

%x = 0:dx:size(Im,1)-1;
%y = 0:dy:size(Im,2)-1;
kx = linspace(-1/dx/2,1/dx/2,size(Im,1));
ky = linspace(-1/dy/2,1/dy/2,size(Im,2));

%
%https://mathworld.wolfram.com/Wiener-KhinchinTheorem.html

%convert to double
Im=double(Im); 

% 2D window
r = 0.1;
w1 = repmat(tukeywin(size(Im,1),r),1,size(Im,2));
w2 = repmat(tukeywin(size(Im,2),r)',size(Im,1),1);
Im = Im.*w1.*w2;

%subtract mean
Im=Im-mean(Im(:)); 

%normalize magnitude
Im=Im/sqrt(sum(Im(:).^2)); 

%compute fft2
I=fft2(Im); 

% auto correlation
% A=real(fftshift(ifft2(I.*conj(I)))); 

%Power
P=abs(fftshift((I.*conj(I))))/numel(Im);
%% calc the radius and theta
rho = zeros(numel(kx),numel(ky));
thetaI = zeros(numel(kx),numel(ky));
for ix= 1 : numel(kx)
    for iy= 1 : numel(ky)
        rho(ix,iy) = sqrt(kx(ix)^2+ky(iy)^2);
        thetaI(ix,iy) = atan2(ky(iy),kx(ix));
    end
end
%% make the interpolation
normk = linspace(0,max(rho(:)),min([numel(kx),numel(ky)]));
out = zeros(size(normk));
for i = 1:numel(normk)
    %out(i) = griddata(rho,thetaI,poissrnd,normk(i),theta,'cubic');
    out(i) = griddata(rho,thetaI,P,normk(i),theta,'nearest'); %#ok<GRIDD> 
end
% return

%Pinter = scatteredInterpolant(rho,thetaI,PI);
%out = Pinter(normk,theta*ones(1,numel(normk)));

%rho = reshape(rho,1,numel(rho))';
%thetaI = reshape(rho,1,numel(thetaI))';
%PI = reshape(P,1,numel(P))';