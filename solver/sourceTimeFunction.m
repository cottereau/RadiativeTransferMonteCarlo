function wave = sourceTimeFunction(kind, parameters)
%SOURCETIMEFUNCTION Physical waveforms, positive-sign Fourier transforms.
% wave = sourceTimeFunction('ricker',struct('f0',10,'delay',0))
% wave = sourceTimeFunction('ormsby',struct('corners',[2 4 8 10],'delay',0))
% wave = sourceTimeFunction('toneburst',struct('fc',10,'cycles',3,'delay',0))
% wave = sourceTimeFunction('sampled',struct('time',t,'values',s))
% Frequencies in Hz in parameters; spectrum argument is rad/s (either sign).
% value(t), derivative(t,order), spectrum(omega) preserve argument shape.
% derivativeEnergy(m) = integral |s^(m)|^2 dt, m=1,2.
% eventTime(m) is the centroid of this derivative power.
% Sampled data define a clamped cubic spline with zero value/slope at both
% ends, extended by zero. Values at both ends must be exactly zero.
% No implicit window, amplitude renormalization, FFT, or toolbox is used.
if nargin < 2, parameters = struct; end
if ~isstruct(parameters) || ~isscalar(parameters)
    error('Fonte:InvalidWaveform','Waveform parameters must be a scalar structure.');
end
kind = validatestring(kind,{'ricker','ormsby','toneburst','sampled'});
switch kind
    case 'ricker', allowed = {'f0','delay'};
    case 'ormsby', allowed = {'corners','delay'};
    case 'toneburst', allowed = {'fc','cycles','delay'};
    case 'sampled', allowed = {'time','values','delay'};
end
if any(~ismember(fieldnames(parameters),allowed))
    error('Fonte:InvalidWaveform','Unknown waveform parameter; set physical amplitude in the source mechanism.');
end
if ~isfield(parameters,'delay'), parameters.delay = 0; end
validateattributes(parameters.delay,{'numeric'},{'real','finite','scalar'});
delay = double(parameters.delay);
wave = struct('kind',kind,'parameters',parameters,'fourierSign',1, ...
    'derivativeEnergy',zeros(1,2),'eventTime',zeros(1,2), ...
    'supportOmega',[0 Inf],'breakpointsOmega',[], ...
    'omegaScale',1,'duration',0);
switch kind
    case 'ricker'
        f = positiveParameter(parameters,'f0'); a = pi*f; omega0 = 2*a;
        wave.value = @(t) rickerDerivative(t,0,a,delay);
        wave.derivative = @(t,m) rickerDerivative(t,m,a,delay);
        wave.spectrum = @(omega) rickerSpectrum(omega,omega0,delay);
        % Exact Gaussian spectral moments, using positive-frequency Parseval.
        for m = 1:2
            wave.derivativeEnergy(m) = ...
                8*omega0^(2*m-1)*gamma(m+2.5)/2^(m+2.5);
        end
        wave.eventTime(:) = delay;
        wave.omegaScale = omega0;
    case 'ormsby'
        if ~isfield(parameters,'corners'), error('Fonte:InvalidWaveform','Missing corners.'); end
        validateattributes(parameters.corners,{'numeric'},{'real','finite','vector','numel',4,'positive'});
        f = double(parameters.corners(:)');
        if any(diff(f)<=0), error('Fonte:InvalidWaveform','Ormsby corners must increase strictly.'); end
        wave.value = @(t) ormsbyDerivative(t,0,f,delay);
        wave.derivative = @(t,m) ormsbyDerivative(t,m,f,delay);
        wave.spectrum = @(omega) ormsbySpectrum(omega,f,delay);
        wave.supportOmega = 2*pi*[f(1) f(4)];
        wave.breakpointsOmega = 2*pi*f;
        wave.omegaScale = 2*pi*f(4);
        [omega,w] = sourceGaussRule(2*pi*f);
        for m = 1:2
            wave.derivativeEnergy(m) = sum(w.*omega.^(2*m).*abs(wave.spectrum(omega)).^2)/pi;
        end
        wave.eventTime(:) = delay;
    case 'toneburst'
        fc = positiveParameter(parameters,'fc');
        cycles = positiveParameter(parameters,'cycles');
        validateattributes(cycles,{'numeric'},{'integer'});
        T = cycles/fc; wc = 2*pi*fc;
        wave.value = @(t) toneDerivative(t,0,wc,cycles,delay);
        wave.derivative = @(t,m) toneDerivative(t,m,wc,cycles,delay);
        wave.spectrum = @(omega) toneSpectrum(omega,wc,cycles,delay);
        % Orthogonality of integer harmonics over the burst duration.
        frequencies = wc*[1,1+1/cycles,1-1/cycles];
        amplitudes = [.5,-.25,-.25];
        for m = 1:2
            wave.derivativeEnergy(m) = T/2*sum(amplitudes.^2.*frequencies.^(2*m));
        end
        wave.eventTime(:) = delay+T/2;
        wave.omegaScale = wc;
        wave.duration = T;
    case 'sampled'
        if ~isfield(parameters,'time') || ~isfield(parameters,'values')
            error('Fonte:InvalidWaveform','Sampled waveforms require time and values.');
        end
        validateattributes(parameters.time,{'numeric'},{'real','finite','vector'});
        validateattributes(parameters.values,{'numeric'},{'real','finite','vector'});
        t = double(parameters.time(:)'); y = double(parameters.values(:)');
        if numel(t)~=numel(y) || numel(t)<3 || any(diff(t)<=0)
            error('Fonte:InvalidWaveform','Use at least three strictly increasing times with matching values.');
        end
        if y(1)~=0 || y(end)~=0 || ~any(y)
            error('Fonte:InvalidWaveform','Sampled pulses need zero endpoints and nonzero interior values.');
        end
        if delay~=0
            error('Fonte:InvalidWaveform','Sampled times are absolute; shift the time vector instead of using delay.');
        end
        % Work relative to the first sample for accurate polynomial evaluation.
        pp = spline(t-t(1),[0 y 0]);
        [x,w] = sourceGaussRule(t-t(1));
        for m = 1:2
            power = splineDerivative(x,m,pp).^2;
            wave.derivativeEnergy(m) = sum(w.*power);
            wave.eventTime(m) = t(1)+sum(w.*x.*power)/wave.derivativeEnergy(m);
        end
        wave.value = @(time) splineDerivative(time-t(1),0,pp);
        wave.derivative = @(time,m) splineDerivative(time-t(1),m,pp);
        wave.spectrum = @(omega) splineSpectrum(omega,pp,t(1));
        wave.omegaScale = pi/min(diff(t));
        wave.duration = t(end)-t(1);
end
if any(~isfinite(wave.derivativeEnergy)) || any(wave.derivativeEnergy<=0)
    error('Fonte:InvalidWaveformEnergy','Waveform derivative energies must be finite and positive.');
end
end

function value = positiveParameter(p,name)
if ~isfield(p,name), error('Fonte:InvalidWaveform','Missing %s.',name); end
value = p.(name);
validateattributes(value,{'numeric'},{'real','finite','scalar','positive'});
value = double(value);
end

function checkOrder(m)
validateattributes(m,{'numeric'},{'real','scalar','integer','>=',0,'<=',2});
end

function y = rickerDerivative(t,m,a,delay)
checkOrder(m); x = a*(t-delay);
switch m
    case 0, y = (1-2*x.^2).*exp(-x.^2);
    case 1, y = a*(4*x.^3-6*x).*exp(-x.^2);
    case 2, y = a^2*(-8*x.^4+24*x.^2-6).*exp(-x.^2);
end
end

function S = rickerSpectrum(omega,w0,delay)
S = (4*sqrt(pi)/w0)*(omega/w0).^2.*exp(-(omega/w0).^2).*exp(1i*omega*delay);
end

function S = ormsbySpectrum(omega,f,delay)
x = abs(omega)/(2*pi);
trapezoid = max(0,min(min((x-f(1))/(f(2)-f(1)),1),(f(4)-x)/(f(4)-f(3))));
S = trapezoid/(f(3)+f(4)-f(1)-f(2)).*exp(1i*omega*delay);
end

function y = ormsbyDerivative(t,m,f,delay)
checkOrder(m); shape = size(t); theta = t(:)-delay;
% Inverse transform of the piecewise linear positive-frequency spectrum.
% Exact exponential moments avoid sinc-derivative cancellation at t=delay.
v = [0 1 1 0]/(f(3)+f(4)-f(1)-f(2));
edges = 2*pi*f; y = zeros(size(theta));
for j = 1:3
    lo = edges(j); h = edges(j+1)-lo;
    % Ascending powers of u: (lo+h*u)^m * (v_j + delta_v*u).
    polynomial = 1;
    for k = 1:m, polynomial = conv(polynomial,[lo h]); end
    polynomial = conv(polynomial,[v(j),v(j+1)-v(j)]);
    moments = sourceExpMoments(-h*theta,numel(polynomial)-1);
    y = y+real((-1i)^m*h*exp(-1i*lo*theta).*(moments*polynomial(:)))/pi;
end
y = reshape(y,shape);
end

function y = toneDerivative(t,m,wc,N,delay)
checkOrder(m); theta = t-delay; T = 2*pi*N/wc;
y = zeros(size(t)); active = theta>=0 & theta<=T;
frequencies = wc*[1,1+1/N,1-1/N]; amplitudes = [.5,-.25,-.25];
for j = 1:3
    y(active) = y(active)+amplitudes(j)*frequencies(j)^m* ...
        sin(frequencies(j)*theta(active)+m*pi/2);
end
% Exact zero endpoints (the analytic limits of all three orders).
y(theta==0 | theta==T) = 0;
end

function S = toneSpectrum(omega,wc,N,delay)
T = 2*pi*N/wc; S = zeros(size(omega));
frequencies = wc*[1,1+1/N,1-1/N]; amplitudes = [.5,-.25,-.25];
for j = 1:3
    S = S+amplitudes(j)/(2i)*(intervalTransform(omega+frequencies(j),T) ...
        -intervalTransform(omega-frequencies(j),T));
end
S = S.*exp(1i*omega*delay);
S(omega==0) = 0;
end

function R = intervalTransform(omega,T)
x = omega*T/2; s = ones(size(x)); ind = x~=0;
s(ind) = sin(x(ind))./x(ind);
R = T*exp(1i*x).*s;
end

function y = splineDerivative(t,m,pp)
checkOrder(m); c = pp.coefs;
for j = 1:m, c = c(:,1:end-1).*(size(c,2)-1:-1:1); end
dpp = mkpp(pp.breaks,c); y = ppval(dpp,t);
y(t<pp.breaks(1) | t>pp.breaks(end)) = 0;
end

function S = splineSpectrum(omega,pp,origin)
shape = size(omega); omega = omega(:); S = zeros(size(omega));
for j = 1:pp.pieces
    h = pp.breaks(j+1)-pp.breaks(j);
    polynomial = fliplr(pp.coefs(j,:)).*h.^(0:3);
    % Block frequency evaluations to bound temporary memory for user grids.
    for first = 1:256:numel(omega)
        ids = first:min(first+255,numel(omega));
        moments = sourceExpMoments(h*omega(ids),3);
        S(ids) = S(ids)+h*exp(1i*omega(ids)*(origin+pp.breaks(j))).*(moments*polynomial(:));
    end
end
S = reshape(S,shape);
end
