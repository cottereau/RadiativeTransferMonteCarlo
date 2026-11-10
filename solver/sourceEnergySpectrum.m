function spectrum = sourceEnergySpectrum(source,wave)
%SOURCEENERGYSPECTRUM Absolute 3D point-source energy, per angular frequency.
% spectrum.P(omega), .S(omega), .total(omega): J/(rad/s), omega>=0.
% .perHzP(f), .perHzS(f): J/Hz. .energy=[EP ES] in J; .eventTime in s.
% Uses the physical mixed source even if source.launchMode is conditional.
% The point-source limit is used; lambda is numerical position jitter only.
% Acoustic sources (sourceAcoustic) return their single mode in .P, .S = 0.
if ~isstruct(source) || ~isscalar(source) || ~isfield(source,'type') || ...
        ~any(strcmp(source.type,{'point_force','moment_tensor','acoustic'}))
    error('Fonte:InvalidSource','Physical energy requires a force, moment or acoustic source.');
end
source.launchMode = 'mixed';
src = physicalSourceForContext(source,3,strcmp(source.type,'acoustic'));
required = {'spectrum','derivativeEnergy','eventTime','omegaScale', ...
    'supportOmega','breakpointsOmega','duration','kind'};
if ~isstruct(wave) || ~isscalar(wave) || ~all(isfield(wave,required)) || ...
        ~isa(wave.spectrum,'function_handle')
    error('Fonte:InvalidWaveform','Use a waveform returned by sourceTimeFunction.');
end
validateattributes(wave.derivativeEnergy,{'numeric'},{'real','finite','positive','numel',2});
if strcmp(src.type,'point_force')
    m = 1; scale = max(abs(src.F));
    strength = sum((src.F/scale).^2);
    logP = 2*log(scale)+log(strength)-log(12*pi^2)-log(src.rho)-3*log(src.vp);
    logS = 2*log(scale)+log(strength)-log(6*pi^2)-log(src.rho)-3*log(src.vs);
elseif strcmp(src.type,'acoustic')
    % E' = omega^2 |S|^2 [ (rho c Q)^2 + |F|^2/3 ] / (4 pi^2 rho c^3);
    % the monopole-dipole cross term integrates to zero over directions.
    m = 1; a = src.rho*src.c*src.Q; scale = max(abs([a src.F]));
    strength = (a/scale)^2+sum((src.F/scale).^2)/3;
    logP = 2*log(scale)+log(strength)-log(4*pi^2)-log(src.rho)-3*log(src.c);
    logS = -Inf;
else
    m = 2; scale = max(abs(src.M(:))); Q = src.M/scale;
    D = sourceMomentDeviator(src.M,scale);
    CP = trace(Q)^2+2*sum(Q(:).^2); CS = 3*sum(D(:).^2);
    logP = 2*log(scale)+log(CP)-log(60*pi^2)-log(src.rho)-5*log(src.vp);
    logS = 2*log(scale)+log(CS)-log(60*pi^2)-log(src.rho)-5*log(src.vs);
end
logNorm = log(pi)+log(wave.derivativeEnergy(m));
energy = exp([logP logS]+logNorm);
if any(~isfinite(energy)) || sum(energy)==0 || ~isfinite(sum(energy))
    error('Fonte:UnrepresentableEnergy','Physical energy is outside the floating-point range; check source units/amplitude.');
end
spectrum = struct('source',src,'waveform',wave,'derivativeOrder',m, ...
    'energy',energy,'totalEnergy',sum(energy),'eventTime',wave.eventTime(m), ...
    'spatialModel','point','numericalJitter',src.lambda, ...
    'densityUnits','J/(rad/s)','energyUnits','J');
spectrum.P = @(omega) evaluate(omega,logP,m,wave);
spectrum.S = @(omega) evaluate(omega,logS,m,wave);
spectrum.total = @(omega) evaluate(omega,logP,m,wave)+evaluate(omega,logS,m,wave);
spectrum.perHzP = @(f) 2*pi*evaluate(2*pi*f,logP,m,wave);
spectrum.perHzS = @(f) 2*pi*evaluate(2*pi*f,logS,m,wave);
spectrum.probabilityDensityOmega = @(omega) evaluate(omega,-logNorm,m,wave);
end

function density = evaluate(omega,logCoefficient,m,wave)
validateattributes(omega,{'numeric'},{'real','finite','nonnegative'});
omega = double(omega); amplitude = abs(wave.spectrum(omega));
density = zeros(size(omega)); active = omega>0 & amplitude>0;
if isinf(logCoefficient) && logCoefficient<0, return; end
density(active) = exp(logCoefficient+2*m*log(omega(active))+2*log(amplitude(active)));
if any(~isfinite(density(:)))
    error('Fonte:UnrepresentableEnergy','Spectral density exceeds the floating-point range.');
end
end
