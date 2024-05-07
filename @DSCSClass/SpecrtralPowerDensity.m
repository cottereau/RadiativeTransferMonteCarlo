ccc

Im=imread('rice.png');

dx = 1;
dy = 1;
theta = [0 pi pi/4 pi/2 3*pi/4 -pi/4 -pi/2 -3*pi/4];
theta = [0 pi];
figure
hold all
for it = 1:numel(theta)
    [P(it,:),k(it,:)] = Image2PSD(Im,dx,dy,theta(it));
    it
    plot(k(it,:),P(it,:),'DisplayName',['\theta ', num2str(theta(it)*pi/180)])
end

%%
figure
hold all
for it = 1:numel(theta)
    plot(k(it,:),P(it,:)/max(P(it,:)),'DisplayName',['\theta ', num2str(theta(it))])
end

return
%%
figure
colormap jet
subplot(1,2,1)
surf(Im)
title('Image')
xlabel('x')
ylabel('y')
axis tight
view(2)
shading flat

subplot(1,2,2)
plot(k,P)
title(['Power Spectral Density @ angle ',num2str(theta)])
xlabel('Wavenumber |K| [1/m]')
ylabel('Power Spectral Density []')
ylim([0 1.1*max(P)])
box on
grid on
return

x = 0:dx:size(Im,1)-1;
y = 0:dy:size(Im,2)-1;
kx = linspace(-1/dx/2,1/dx/2,size(Im,1));
ky = linspace(-1/dy/2,1/dy/2,size(Im,2));

figure
colormap jet
subplot(2,2,1)
surf(x,y,Im)
title('Image')
xlabel('x')
ylabel('y')
axis tight
view(2)
shading flat

%convert to double
Im=double(Im); 
%subtract mean
Im=Im-mean(Im(:)); 
%normalize magnitude
Im=Im/sqrt(sum(Im(:).^2)); 
%compute fft2
I=fftshift(fft2(Im)); 

subplot(2,2,3)
colormap jet
surf(kx,ky,abs(I))
title('Image abs fft')
xlabel('k_x')
ylabel('k_y')
axis tight
caxis([0 5])
view(2)
shading flat
subplot(2,2,4)
surf(kx,ky,angle(I))
title('Image angle fft')
colormap jet
xlabel('k_x')
ylabel('k_y')
axis tight
view(2)
shading flat

% autocorrelation
A=real(ifft2(I.*conj(I))); 
%Power
P=(I.*conj(I))/numel(Im);

subplot(2,2,2)
surf(kx,ky,P)
title('PSD')
xlabel('k_x')
ylabel('k_y')
axis tight
view(2)
shading flat
caxis([0 0.0005])
colormap jet

%% calc the radius and theta

for ix= 1 : numel(kx)
    for iy= 1 : numel(ky)
        r(ix,iy) = sqrt(kx(ix)^2+ky(iy)^2);
        theta(ix,iy) = atan2(ky(iy),kx(ix));
    end
end
