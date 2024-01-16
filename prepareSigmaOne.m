function [Sigma,Sigma_prime,invcdf] = prepareSigmaOne(sigma,d)
Nth = 1000;
xth = linspace(0,pi,Nth);
if d==2
    Sigma = 2*integral(sigma,0,pi);
    Sigma_prime = 2*integral(@(th)sigma(th).*cos(th),0,pi);
    sigmaNorm = @(th) (2/Sigma)*sigma(th);
elseif d==3
    Sigma = 2*pi*integral(@(th)sigma(th).*sin(th),0,pi);
    Sigma_prime = 2*pi*integral(@(th)sigma(th).*sin(th).*cos(th),0,pi);
    sigmaNorm = @(th) (2*pi/Sigma)*sin(th).*sigma(th);
end
pdf = sigmaNorm(xth);
cdf = cumsum(pdf)*(pi/Nth);
ind = [find(diff(cdf)>0) Nth];
ind(1) = 1;
invcdf = griddedInterpolant(cdf(ind),xth(ind));
end
