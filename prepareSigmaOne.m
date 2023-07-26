function [Sigma,invcdf] = prepareSigmaOne(sigma,d)
Nth = 1000;
% Question : why the dependence of sigma/Sigma on q/k is not considered?
% lc (correlation distance) and the correlation model should also be given
% as input ...
% sigma is defined rather a function of cos(\thet)=\xi instead of the scattering angle theta itself. 
% Hence, the intergation will be done over (-1,+1).
if d==2
    Sigma = 2*integral(sigma,0,pi);
    sigmaNorm = @(th) 2*sigma(th)/Sigma;
    ind = 1:Nth;
elseif d==3
    Sigma = 2*pi*integral(@(th)sigma(th).*sin(th),0,pi);
    sigmaNorm = @(th) sigma(th).*sin(th)/Sigma;
    ind = [1 3:Nth-2 Nth];
end
xth = linspace(0,pi,Nth);
pdf = sigmaNorm(xth);
cdf = cumsum(pdf)*(pi/Nth);
invcdf = griddedInterpolant(cdf(ind),xth(ind));
end
