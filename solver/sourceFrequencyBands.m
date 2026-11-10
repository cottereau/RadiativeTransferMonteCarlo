function bands = sourceFrequencyBands(spectrum,opts)
%SOURCEFREQUENCYBANDS Positive quadrature with bounded spectral refinement.
% bands.omega and weightOmega are rad/s; frequencyHz is for MaterialClass.
% energyP/S are physical weighted node energies; never renormalized.
% opts: tailTolerance (1e-3), quadratureTolerance (1e-6), maxNodes (256),
% maxTailExpansions (6). Errors on budget exhaustion instead of large runs.
% Convergence here concerns the source spectrum ONLY, not receiver responses.
if nargin < 2, opts = struct; end
if ~isstruct(opts) || ~isscalar(opts), error('Fonte:InvalidBands','Options must be a scalar structure.'); end
defaults = struct('tailTolerance',1e-3,'quadratureTolerance',1e-6, ...
    'maxNodes',256,'maxTailExpansions',6);
names = fieldnames(defaults);
if any(~ismember(fieldnames(opts),names)), error('Fonte:InvalidBands','Unknown frequency-band option.'); end
for j = 1:numel(names)
    name = names{j};
    if ~isfield(opts,name), opts.(name) = defaults.(name); end
end
validateattributes(opts.tailTolerance,{'numeric'},{'real','finite','scalar','>',0,'<',1});
validateattributes(opts.quadratureTolerance,{'numeric'},{'real','finite','scalar','>',0,'<',1});
validateattributes(opts.maxNodes,{'numeric'},{'real','finite','scalar','integer','>=',8});
validateattributes(opts.maxTailExpansions,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
if ~isstruct(spectrum) || ~isfield(spectrum,'probabilityDensityOmega') || ~isfield(spectrum,'waveform')
    error('Fonte:InvalidBands','Use a spectrum returned by sourceEnergySpectrum.');
end
wave = spectrum.waveform; scale = wave.omegaScale;
compact = isfinite(wave.supportOmega(2));
if compact
    bounds = wave.supportOmega/scale;
else
    bounds = [0 2];
    if strcmp(wave.kind,'ricker'), bounds(2) = 4; end
end
for expansion = 0:opts.maxTailExpansions
    breaks = unique([bounds(1),wave.breakpointsOmega/scale,bounds(2)]);
    breaks = breaks(breaks>=bounds(1) & breaks<=bounds(2));
    % Resolve finite-duration spectral oscillations before error comparison.
    panels = max(1,ceil(diff(breaks)*scale*wave.duration/pi));
    if isscalar(panels), panels = max(panels,4); end
    [~,~,retained] = integratePanels(panels);
    while true
        if 8*sum(panels)>opts.maxNodes
            error('Fonte:BandBudget','Source quadrature needs more than maxNodes=%d; relax tolerance or explicitly raise the budget.',opts.maxNodes);
        end
        panels = 2*panels;
        [x,w,fine] = integratePanels(panels);
        quadratureError = abs(fine-retained);
        retained = fine;
        if quadratureError<=opts.quadratureTolerance, break; end
    end
    if retained>1+max(4*quadratureError,1e-10)
        error('Fonte:SpectrumNormalization','Integrated spectrum exceeds the independent total-energy normalization.');
    end
    omitted = max(0,1-retained);
    if compact
        if abs(1-retained)>max(4*quadratureError,opts.quadratureTolerance)
            error('Fonte:SpectrumNormalization','Compact spectrum disagrees with its total energy.');
        end
        break
    end
    if omitted+quadratureError<=opts.tailTolerance, break; end
    if expansion==opts.maxTailExpansions
        error('Fonte:TailBudget','Tail tolerance was not reached within maxTailExpansions.');
    end
    bounds(2) = 2*bounds(2);
end
omega = scale*x; weights = scale*w;
energyP = weights.*spectrum.P(omega); energyS = weights.*spectrum.S(omega);
retainedEnergy = [sum(energyP),sum(energyS)];
if compact, omitted = 0; end % exact spectral support; quadrature error separate
bands = struct('omega',omega,'frequencyHz',omega/(2*pi),'weightOmega',weights, ...
    'energyP',energyP,'energyS',energyS,'totalEnergy',spectrum.totalEnergy, ...
    'modalEnergy',spectrum.energy,'retainedEnergy',retainedEnergy, ...
    'retainedTotalEnergy',sum(retainedEnergy),'omittedEnergy',omitted*spectrum.energy, ...
    'omittedTotalEnergy',omitted*spectrum.totalEnergy, ...
    'omittedFraction',omitted,'quadratureErrorEstimate',quadratureError, ...
    'normalizationResidual',1-retained,'omegaRange',scale*bounds, ...
    'eventTime',spectrum.eventTime,'nodeCount',numel(omega), ...
    'tailExpansions',expansion,'options',opts, ...
    'responseConvergence','not_checked','energyUnits','J');

    function [nodes,weights,fraction] = integratePanels(counts)
        if 4*sum(counts)>opts.maxNodes
            error('Fonte:BandBudget','Resolving this waveform requires more than maxNodes=%d.',opts.maxNodes);
        end
        edges = breaks(1);
        for k = 1:numel(counts)
            local = linspace(breaks(k),breaks(k+1),counts(k)+1);
            edges = [edges local(2:end)]; %#ok<AGROW>
        end
        [nodes,weights] = sourceGaussRule(edges);
        fraction = sum(weights.*(scale*spectrum.probabilityDensityOmega(scale*nodes)));
    end
end
