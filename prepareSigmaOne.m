function [Sigma,invcdf] = prepareSigmaOne(sigma,d)
Nth = 1000;
% Question : why the dependence of sigma/Sigma on q/k is not considered?
% lc (correlation distance) and the correlation model should also be given
% as input ...
% sigma is defined rather a function of cos(\thet)=\xi instead of the scattering angle theta itself. 
% Hence, the intergation will be done over (-1,+1).
xth = linspace(0,pi,Nth);
if any(sigma(xth)<0)
    warning('sigma function should be positive. This does not seem to be the case')
    sigma = @(th) max(0,sigma(th));
end
if d==2
    Sigma = 2*integral(sigma,0,pi);
    sigmaNorm = @(th) (2/Sigma)*sigma(th);
elseif d==3
    Sigma = 2*pi*integral(@(th)sigma(th).*sin(th),0,pi);
    sigmaNorm = @(th) (2*pi/Sigma)*sin(th).*sigma(th);
end
pdf = sigmaNorm(xth);
cdf = cumsum(pdf)*(pi/Nth);
ind = [find(diff(cdf)>0) Nth];
ind(1) = 1;
invcdf = griddedInterpolant(cdf(ind),xth(ind));
end
