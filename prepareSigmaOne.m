function [Sigma,invcdf] = prepareSigmaOne(sigma)
Nth = 1000;
% Question : why the dependence of sigma/Sigma on q/k is not considered?
% lc (correlation distance) and the correlation model should also be given
% as input ...
% sigma is defined rather a function of cos(\thet)=\xi instead of the scattering angle theta itself. 
% Hence, the intergation will be done over (-1,+1).
Sigma = integral(sigma,0,2*pi);
sigmaNorm = @(th) sigma(th)/Sigma;
xth = linspace(0,2*pi,Nth);
pdf = sigmaNorm(xth);
cdf = cumsum(pdf)*2*pi/Nth;
invcdf = griddedInterpolant(cdf,xth);
end
