function TPCF = PolycrystalTPCF(modelName,modelParameter)
%POLYCRYSTALTPCF Construct a dimensional spectral TPCF for a 3D polycrystal.
%
% The returned TPCF.spectrum(q) follows the Fourier convention
%
%   etaTilde(q) = (2*pi)^(-3) integral eta(r) exp(-i*q.r) dr
%
% and consequently has units of length cubed. Available models are:
%
%   'exponential'    eta(r) = exp(-r/a), where modelParameter is a.
%   'spherical'      Equal spherical grains, where modelParameter is D.
%   'Arguelles'      Exponential kernel averaged over the GSD.
%   'Sha'            Spherical kernel averaged over the GSD.
%   'ShengKhazaie'   Volume-weighted spherical kernel averaged over the GSD.
%
% For the GSD-dependent models, modelParameter is a scalar structure. The
% supported distribution types are 'lognormal', 'gamma', 'weibull', and
% 'truncatedNormal':
%
%   grainSizeDistribution = struct( ...
%       'type','gamma', ...
%       'meanDiameter',meanDiameter, ...
%       'standardDeviation',standardDeviation);
%
% The mean and standard deviation always describe the grain diameter 
% distribution. In particular, 'truncatedNormal' denotes an underlying 
% Gaussian distribution conditioned on D > 0; the parameters of that 
% underlying Gaussian are calculated internally.

if ~(ischar(modelName) || (isstring(modelName) && isscalar(modelName)))
    error('MaterialClass:InvalidPolycrystalTPCFName', ...
          'The polycrystal TPCF model must be a character vector or scalar string.');
end

key = lower(strtrim(char(modelName)));
key = regexprep(key,'[\s_-]',''); % Remove spaces, underscores, and hyphens

switch key
    case {'exp','exponential'}
        correlationLength = modelParameter;
        TPCF = struct( ...
            'model','exponential', ...
            'parameters',struct('correlationLength',correlationLength), ...
            'spectrum',@(q) correlationLength^3./(pi^2*(1+(correlationLength.*q).^2).^2));

    case {'spherical','sphericalmonodisperse','monodispersespherical'}
        grainDiameter = modelParameter;
        TPCF = struct( ...
            'model','spherical', ...
            'parameters',struct('grainDiameter',grainDiameter), ...
            'spectrum',@(q) grainDiameter^3.*sphericalSpectrum(q.*grainDiameter));

    case {'arguelles','arguellesturner'}
        [parameters,diameters,probabilityWeights,thirdMoment] = grainDiameterQuadrature(modelParameter,3);
        spectralWeights = thirdMoment.*probabilityWeights;
        parameters.thirdMoment = thirdMoment;
        TPCF = struct( ...
            'model','Arguelles', ...
            'parameters',parameters, ...
            'spectrum',@(q) averagedSpectrum( ...
                q,diameters,spectralWeights,'exponential'));

    case 'sha'
        [parameters,diameters,probabilityWeights,thirdMoment] = grainDiameterQuadrature(modelParameter,3);
        spectralWeights = thirdMoment.*probabilityWeights;
        parameters.thirdMoment = thirdMoment;
        TPCF = struct( ...
            'model','Sha', ...
            'parameters',parameters, ...
            'spectrum',@(q) averagedSpectrum(q,diameters,spectralWeights,'spherical'));

    case {'shengkhazaie','sk'}
        [parameters,diameters,probabilityWeights,sixthMoment] = grainDiameterQuadrature(modelParameter,6);
        thirdMoment = grainDiameterMoment(parameters,3);
        spectralWeights = sixthMoment/thirdMoment.*probabilityWeights;
        parameters.thirdMoment = thirdMoment;
        parameters.sixthMoment = sixthMoment;
        TPCF = struct( ...
            'model','ShengKhazaie', ...
            'parameters',parameters, ...
            'spectrum',@(q) averagedSpectrum(q,diameters,spectralWeights,'spherical'));

    otherwise
        error('MaterialClass:UnknownPolycrystalTPCF', ...
             ['Unknown polycrystal TPCF model "%s". Available models are ', ...
              '''exponential'', ''spherical'', ''Arguelles'', ''Sha'', ', ...
              'and ''ShengKhazaie''.'],char(modelName));
end
end

function [parameters,diameters,probabilityWeights,diameterMoment] = ...
        grainDiameterQuadrature(grainSizeDistribution,diameterPower)
%GRAINDIAMETERQUADRATURE Construct a quadrature for a moment-weighted GSD.
%
% The spectral TPCF models require diameter averages of the form
%
%   integral D^p*K(q*D)*p_D(D) dD,
%
% where p_D is the grain-diameter probability density, p is
% diameterPower, and K is the model-dependent spectral kernel. 
% Introducing the normalized moment-weighted density
%
%   p_p(D) = D^p*p_D(D)/E[D^p]
%
% gives the exact change-of-measure identity
%
%   E[D^p] * integral K(q*D)*p_p(D) dD
%        = integral D^p*K(q*D)*p_D(D) dD.
%
% Only the remaining integral with respect to p_p is approximated:
%
%   diameterMoment * sum_j probabilityWeights(j)*K(q*diameters(j)).
%
% This is not the generally false factorization
% E[D^p*K(q*D)] = E[D^p]*E[K(q*D)]. The expectation of K on the right-hand
% side is taken with respect to the new density p_p, not the original GSD.
%
% Thus, diameters are deterministic quadrature nodes, not random grain
% samples. Likewise, probabilityWeights are normalized quadrature weights
% for p_p(D), not probabilities assigned to bins of the original GSD. They
% sum to one. This moment weighting places more quadrature resolution on
% the larger grains that contribute most strongly to the spectrum.
%
% No Monte Carlo sampling is used. The deterministic integration rule
% depends on the selected GSD:
%   lognormal      - composite trapezoidal rule after transforming the
%                    log-diameter to a standard-normal variable on [-10,10]
%   gamma          - generalized Gauss-Laguerre quadrature, constructed
%                    with the Golub-Welsch eigenvalue method
%   weibull        - composite trapezoidal rule in scaled diameter, with
%                    the upper limit obtained from the moment-weighted tail
%   truncatedNormal - composite trapezoidal rule in diameter, or in a
%                    rescaled tail variable when the truncation is severe
%   monodisperse   - exact one-node rule
%
% Inputs
%   grainSizeDistribution : scalar structure with the fields
%       type              - 'lognormal', 'gamma', 'weibull', or
%                           'truncatedNormal'
%       meanDiameter      - mean grain diameter of the positive GSD [m]
%       standardDeviation - standard deviation of grain diameter [m]
%   diameterPower : exponent p in the required diameter moment.
%                   p = 3 for the Arguelles and Sha models.
%                   p = 6 for the numerator of the ShengKhazaie model.
%
% Outputs
%   parameters         : canonical GSD name, prescribed statistics,
%                        distribution-specific parameters, and quadrature
%                        order
%   diameters          : row vector of positive quadrature nodes D_j [m]
%   probabilityWeights : row vector of normalized quadrature weights w_j
%                        for D^p*p_D(D)/E[D^p]
%   diameterMoment     : analytical moment E[D^p] [m^p]
%
% When standardDeviation is zero, the GSD is monodisperse and the
% quadrature reduces exactly to one node at meanDiameter with unit weight.

% Check that the GSD container and its required fields are present.
if ~isstruct(grainSizeDistribution) || ~isscalar(grainSizeDistribution)
    error('MaterialClass:InvalidGrainSizeDistribution', ...
         ['A GSD-dependent TPCF requires a scalar structure containing ', ...
          'type, meanDiameter, and standardDeviation.']);
end

requiredFields = {'type','meanDiameter','standardDeviation'};
for fieldIndex = 1:numel(requiredFields)
    if ~isfield(grainSizeDistribution,requiredFields{fieldIndex})
        error('MaterialClass:InvalidGrainSizeDistribution', ...
              'The grain-size distribution is missing the field %s.', ...
            requiredFields{fieldIndex});
    end
end

distributionName = grainSizeDistribution.type;

% Convert accepted spelling variants to one internal distribution name.
distributionKey = lower(strtrim(char(distributionName)));
distributionKey = regexprep(distributionKey,'[\s_-]','');

switch distributionKey
    case 'lognormal'
        distName = 'lognormal';
    case 'gamma'
        distName = 'gamma';
    case 'weibull'
        distName = 'weibull';
    case {'truncatednormal','truncatedgaussian'}
        distName = 'truncatedNormal';
    otherwise
        error('MaterialClass:UnsupportedGrainSizeDistribution', ...
             ['Unknown grain-size distribution "%s". Available ', ...
              'distributions are ''lognormal'', ''gamma'', ''weibull'', ', ...
              'and ''truncatedNormal''.'],char(distributionName));
end

meanDiameter = grainSizeDistribution.meanDiameter;
standardDeviation = grainSizeDistribution.standardDeviation;

% The coefficient of variation controls the breadth of the GSD and the
% quadrature order required to resolve it.
coeffOfVariation = standardDeviation/meanDiameter;

if standardDeviation == 0
    % Exact monodisperse limit: D is always equal to meanDiameter.
    diameters = meanDiameter;
    probabilityWeights = 1;
    quadratureOrder = 1;
    diameterMoment = meanDiameter^diameterPower;
else
    % Construct distribution-specific nodes and normalized weights for the
    % moment-weighted density D^p*p_D(D)/E[D^p].
    quadratureOrder = selectQuadratureOrder(coeffOfVariation,distName);
    switch distName
        case 'lognormal'
            [diameters,probabilityWeights,diameterMoment,parameters] = ...
                lognormalQuadrature(meanDiameter,coeffOfVariation,diameterPower,quadratureOrder);

        case 'gamma'
            [diameters,probabilityWeights,diameterMoment,parameters] = ...
                gammaQuadrature(meanDiameter,coeffOfVariation,diameterPower,quadratureOrder);

        case 'weibull'
            [diameters,probabilityWeights,diameterMoment,parameters] = ...
                weibullQuadrature(meanDiameter,coeffOfVariation,diameterPower,quadratureOrder);

        case 'truncatedNormal'
            if coeffOfVariation >= 1
                error('MaterialClass:InvalidTruncatedNormalGSD', ...
                     ['A normal distribution truncated at D = 0 has a ', ...
                      'coefficient of variation strictly smaller than 1.']);
            end
            [diameters,probabilityWeights,diameterMoment,parameters] = ...
                truncatedNormalQuadrature(meanDiameter,coeffOfVariation,diameterPower,quadratureOrder);
    end
end

% Give all model branches the same row-vector convention.
diameters = reshape(diameters,1,[]);
probabilityWeights = reshape(probabilityWeights,1,[]);

% Store the canonical GSD description for inspection and for evaluating
% other diameter moments without reconstructing the quadrature.
if standardDeviation == 0
    parameters = struct();
end
parameters.distribution = distName;
parameters.meanDiameter = meanDiameter;
parameters.standardDeviation = standardDeviation;
parameters.coeffOfVariation = coeffOfVariation;
parameters.quadratureOrder = quadratureOrder;
end

function diameterMoment = grainDiameterMoment(parameters,diameterPower)
% Evaluate a diameter moment without rebuilding the quadrature.

if parameters.standardDeviation == 0
    diameterMoment = parameters.meanDiameter^diameterPower;
    return
end

switch parameters.distribution
    case 'lognormal'
        diameterMoment = exp(diameterPower*parameters.underlyingLogMean ...
            + 0.5*diameterPower^2*parameters.underlyingLogStandardDeviation^2);
    case 'gamma'
        diameterMoment = exp(diameterPower*log(parameters.scale) ...
            + gammaln(parameters.shape+diameterPower) ...
            - gammaln(parameters.shape));
    case 'weibull'
        diameterMoment = exp(diameterPower*log(parameters.scale) ...
            + gammaln(1+diameterPower/parameters.shape));
    case 'truncatedNormal'
        standardMoments = truncatedNormalStandardMoments( ...
            parameters.standardizedTruncation,diameterPower);
        diameterMoment = ...
            parameters.underlyingNormalStandardDeviation^diameterPower ...
            *standardMoments(diameterPower+1);
end
end

function quadratureOrder = selectQuadratureOrder(coeffOfVariation,distributionName)
% Use more diameter samples when the distribution becomes broader.

if strcmp(distributionName,'gamma')
    if coeffOfVariation <= 0.5
        quadratureOrder = 128;
    elseif coeffOfVariation <= 0.8
        quadratureOrder = 256;
    else
        quadratureOrder = 512;
    end
elseif any(strcmp(distributionName,{'weibull','truncatedNormal'}))
    if coeffOfVariation <= 0.5
        quadratureOrder = 2049;
    else
        quadratureOrder = 4097;
    end
else
    if coeffOfVariation <= 0.5
        quadratureOrder = 513;
    elseif coeffOfVariation <= 0.8
        quadratureOrder = 2049;
    else
        quadratureOrder = 4097;
    end
end
end

function [diameters,weights,diameterMoment,parameters] = ...
         lognormalQuadrature(meanDiameter,coeffOfVariation,diameterPower,quadratureOrder)
% Evaluate the moment-weighted lognormal in its underlying normal variable.

logVariance = log1p(coeffOfVariation^2);
logStandardDeviation = sqrt(logVariance);
logMean = log(meanDiameter) - 0.5*logVariance;

nodes = linspace(-10,10,quadratureOrder);
weights = exp(-0.5*nodes.^2)/sqrt(2*pi);
weights([1 end]) = 0.5*weights([1 end]);
weights = weights/sum(weights);

diameterMoment = exp(diameterPower*logMean + 0.5*diameterPower^2*logVariance);
tiltedLogMean = logMean + diameterPower*logVariance;
diameters = exp(tiltedLogMean + logStandardDeviation*nodes);

parameters = struct( ...
    'underlyingLogMean',logMean, ...
    'underlyingLogStandardDeviation',logStandardDeviation);
end

function [diameters,weights,diameterMoment,parameters] = ...
          gammaQuadrature(meanDiameter,coeffOfVariation,diameterPower,quadratureOrder)
% A D^p-weighted Gamma distribution is another Gamma distribution.

shape = 1/coeffOfVariation^2;
scale = meanDiameter/shape;
diameterMoment = exp(diameterPower*log(scale) + gammaln(shape+diameterPower)-gammaln(shape));

[transformedDiameter,weights] = generalizedLaguerreQuadrature( ...
    shape+diameterPower,quadratureOrder);
diameters = scale*transformedDiameter;

parameters = struct('shape',shape,'scale',scale);
end

function [diameters,weights,diameterMoment,parameters] = ...
          weibullQuadrature(meanDiameter,coeffOfVariation,diameterPower,quadratureOrder)
% Transform the moment-weighted Weibull distribution to a Gamma variable.

targetLogMomentRatio = log1p(coeffOfVariation^2);
shapeEquation = @(shape) gammaln(1+2/shape) - 2*gammaln(1+1/shape) - targetLogMomentRatio;

lowerShape = 0.1;
while shapeEquation(lowerShape) < 0
    lowerShape = lowerShape/2;
end
upperShape = 10;
while shapeEquation(upperShape) > 0
    upperShape = 2*upperShape;
    if upperShape > 1e8
        error('MaterialClass:InvalidWeibullGSD', ...
              'Could not determine the Weibull shape parameter.');
    end
end

shape = fzero(shapeEquation,[lowerShape upperShape]);
scale = meanDiameter/exp(gammaln(1+1/shape));
diameterMoment = exp(diameterPower*log(scale) + gammaln(1+diameterPower/shape));

tiltedGammaShape = 1 + diameterPower/shape;
maximumTransformedDiameter = gammaincinv(1-1e-14,tiltedGammaShape);
scaledDiameter = linspace(0,maximumTransformedDiameter^(1/shape),quadratureOrder);
weights = scaledDiameter.^(shape+diameterPower-1).*exp(-scaledDiameter.^shape);
weights([1 end]) = 0.5*weights([1 end]);
weights = weights/sum(weights);
diameters = scale*scaledDiameter;

parameters = struct('shape',shape,'scale',scale);
end

function [nodes,weights] = generalizedLaguerreQuadrature(gammaShape,quadratureOrder)
% Integrate expectations under a unit-scale Gamma distribution.
%
% The eigenvalues of this symmetric tridiagonal Jacobi matrix are the
% generalized Gauss-Laguerre nodes. Squared first components of its
% normalized eigenvectors give weights already normalized for a Gamma PDF.

index = (1:quadratureOrder).';
diagonal = 2*index-1+(gammaShape-1);
offDiagonalIndex = (1:quadratureOrder-1).';
offDiagonal = sqrt(offDiagonalIndex .* (offDiagonalIndex+gammaShape-1));
jacobiMatrix = diag(diagonal) + diag(offDiagonal,1) + diag(offDiagonal,-1);

[eigenvectors,nodes] = eig(jacobiMatrix,'vector');
[nodes,order] = sort(nodes);
weights = eigenvectors(1,order).^2;
weights = weights/sum(weights);

nodes = reshape(nodes,1,[]);
weights = reshape(weights,1,[]);
end

function [diameters,weights,diameterMoment,parameters] = ...
          truncatedNormalQuadrature(meanDiameter,coeffOfVariation,diameterPower,quadratureOrder)
% Gaussian diameters conditioned on D > 0 with prescribed final moments.

alpha = solveTruncatedNormalAlpha(coeffOfVariation);
standardMoments = truncatedNormalStandardMoments(alpha,diameterPower);
underlyingStandardDeviation = meanDiameter/standardMoments(2);
underlyingMean = -alpha*underlyingStandardDeviation;
diameterMoment = underlyingStandardDeviation^diameterPower * standardMoments(diameterPower+1);

if alpha <= 5
    tailAtTruncation = 0.5*erfc(alpha/sqrt(2));
    remainingTail = 1e-14*tailAtTruncation;
    maximumNormalVariable = sqrt(2)*erfcinv(2*remainingTail);
    maximumDiameter = underlyingStandardDeviation *(maximumNormalVariable-alpha);
    diameters = linspace(0,maximumDiameter,quadratureOrder);
    normalVariable = diameters/underlyingStandardDeviation+alpha;
    baseWeights = exp(-0.5*normalVariable.^2);
    baseWeights([1 end]) = 0.5*baseWeights([1 end]);
else
    % With t=alpha*(z-alpha), the scaled tail density is well conditioned
    % even when the Gaussian truncation point lies far into its right tail.
    transformedVariable = linspace(0,40,quadratureOrder);
    baseWeights = exp(-transformedVariable ...
        - 0.5*(transformedVariable/alpha).^2);
    baseWeights([1 end]) = 0.5*baseWeights([1 end]);
    diameters = underlyingStandardDeviation*transformedVariable/alpha;
end

scaledDiameter = diameters/meanDiameter;
weights = baseWeights.*scaledDiameter.^diameterPower;
weights = weights/sum(weights);

parameters = struct( ...
    'lowerTruncation',0, ...
    'underlyingNormalMean',underlyingMean, ...
    'underlyingNormalStandardDeviation',underlyingStandardDeviation, ...
    'standardizedTruncation',alpha);
end

function alpha = solveTruncatedNormalAlpha(coeffOfVariation)
% Recover the underlying Gaussian from the moments after truncation at zero.

equation = @(value) truncatedNormalCV(value)-coeffOfVariation;
lowerBound = -max(10,2/coeffOfVariation);
upperBound = 5;
while equation(upperBound) < 0
    upperBound = 2*upperBound;
    if upperBound > 1e4
        error('MaterialClass:InvalidTruncatedNormalGSD', ...
              'Could not determine the underlying truncated-normal parameters.');
    end
end
alpha = fzero(equation,[lowerBound upperBound]);
end

function coeffOfVariation = truncatedNormalCV(alpha)
% Coefficient of variation of Z-alpha conditional on Z > alpha.

moments = truncatedNormalStandardMoments(alpha,2);
variance = max(0,moments(3)-moments(2)^2);
coeffOfVariation = sqrt(variance)/moments(2);
end

function moments = truncatedNormalStandardMoments(alpha,maximumOrder)
% Moments of Y=Z-alpha conditional on Z>alpha for standard-normal Z.

moments = zeros(1,maximumOrder+1);
if alpha <= 5
    tailProbability = 0.5*erfc(alpha/sqrt(2));
    normalDensity = exp(-0.5*alpha^2)/sqrt(2*pi);
    excessIntegrals = zeros(1,maximumOrder+1);
    excessIntegrals(1) = tailProbability;
    if maximumOrder >= 1
        excessIntegrals(2) = normalDensity-alpha*tailProbability;
    end
    for order = 2:maximumOrder
        excessIntegrals(order+1) = (order-1)*excessIntegrals(order-1) ...
            - alpha*excessIntegrals(order);
    end
    moments = excessIntegrals/tailProbability;
else
    % The exponential factor exp(alpha^2/2) cancels from every ratio.
    density = @(value) exp(-value-0.5*(value/alpha).^2);
    normalization = integral(density,0,Inf);
    for order = 0:maximumOrder
        scaledMoment = integral( ...
            @(value) value.^order.*density(value),0,Inf);
        moments(order+1) = scaledMoment/normalization/alpha^order;
    end
end
end

function value = averagedSpectrum(q,diameters,spectralWeights,kernel)
% Evaluate all requested wavenumbers in one vectorized diameter average.

if ~isnumeric(q) || ~isreal(q) || any(~isfinite(q),'all')
    error('MaterialClass:InvalidSpectralWavenumber', ...
          'The spectral wavenumber q must contain finite real values.');
end

originalSize = size(q);
dimensionlessWavenumber = abs(q(:))*diameters;

switch kernel
    case 'exponential'
        % Fourier transform of exp(-2*r/D), with z=q*D.
        kernelValue = 1./(8*pi^2) ./ (1 + dimensionlessWavenumber.^2/4).^2;
    case 'spherical'
        kernelValue = sphericalSpectrum(dimensionlessWavenumber);
    otherwise
        error('MaterialClass:UnknownPolycrystalKernel', ...
              'Unknown internal polycrystal spectral kernel.');
end

value = reshape(sum(kernelValue.*spectralWeights,2),originalSize);
end

function value = sphericalSpectrum(z)
% Dimensionless spectrum of the 3D equal-sphere TPCF

value = zeros(size(z));
small = abs(z) < 0.1;

zSmall = z(small);
value(small) = 1/(48*pi^2).*( 1 - zSmall.^2/20 + 3*zSmall.^4/2800);

zLarge = z(~small);
value(~small) = 3.*(2*sin(zLarge/2) - zLarge.*cos(zLarge/2)).^2 ./ (pi^2*zLarge.^6);
end
