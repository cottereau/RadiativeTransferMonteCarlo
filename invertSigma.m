function invcdf = invertSigma(sigma)
Nth = 1000;
xth = linspace(0,2*pi,Nth);
pdf = sigma(xth);
cdf = cumsum(pdf)*2*pi/Nth;
invcdf = griddedInterpolant(cdf,xth);
end
