function out = composeSourceBands(bands,responses)
%COMPOSESOURCEBANDS Sum physical modal band energies times unit responses.
% responses{j,1/2} is the P/S-injected radiativeTransfer observation.
% A response is required exactly when its weighted band energy is positive.
if ~isstruct(bands) || ~isscalar(bands) || ...
        ~all(isfield(bands,{'energyP','energyS','eventTime'}))
    error('Fonte:InvalidBands','Bands need energyP, energyS and eventTime.');
end
p = double(bands.energyP(:)); s = double(bands.energyS(:)); n = numel(p);
if numel(s)~=n || ~iscell(responses) || ~isequal(size(responses),[n 2]) || ...
        any(~isfinite([p;s])) || any([p;s]<0)
    error('Fonte:InvalidBands','Band energies and response dimensions are invalid.');
end
validateattributes(bands.eventTime,{'numeric'},{'real','finite','scalar'});
first = [];
for j=1:n
    for beta=1:2
        e = [p(j) s(j)];
        r = responses{j,beta};
        if e(beta)==0
            if ~isempty(r), error('Fonte:UnexpectedResponse','Zero-energy channels must be empty.'); end
            continue
        end
        if ~isstruct(r) || ~isscalar(r) || ~isfield(r,'unitInjectedEnergy') || ...
                ~isequal(r.unitInjectedEnergy,true) || ~isfield(r,'energyDensity') || ...
                ~isfield(r,'energy') || ~isfield(r,'aliveEnergy') || ...
                ~isfield(r,'cellMeasure') || ~isfield(r,'t') || ...
                ~isfield(r,'eventTime') || ~isfield(r,'modeLabels')
            error('Fonte:InvalidResponse','Each active channel needs a unit-energy observation.');
        end
        if isempty(first), first = r; end
    end
end
if isempty(first), error('Fonte:ZeroBandEnergy','No positive-energy channel can be composed.'); end
acoustic = isequal(first.modeLabels,{'A'}) && all(s==0);
if ~(isequal(first.modeLabels,{'P','S'}) || acoustic) || first.eventTime~=bands.eventTime
    error('Fonte:ResponseMismatch','Expected elastic P/S (or acoustic) modes at the band event time.');
end
out = first;
out.energyDensity = zeros(size(first.energyDensity));
out.energy = zeros(size(first.energy));
out.aliveEnergy = zeros(size(first.aliveEnergy));
for j=1:n
    for beta=1:2
        e = [p(j) s(j)]; if e(beta)==0, continue; end
        r = responses{j,beta};
        if ~isequal(size(r.energyDensity),size(first.energyDensity)) || ...
                ~isequal(size(r.energy),size(first.energy)) || ...
                ~isequal(size(r.aliveEnergy),size(first.aliveEnergy)) || ...
                ~isequal(r.cellMeasure,first.cellMeasure) || ...
                ~isequal(r.t,first.t) || ~isequal(r.modeLabels,first.modeLabels) || ...
                r.eventTime~=bands.eventTime
            error('Fonte:ResponseMismatch','Unit responses use different grids or event times.');
        end
        out.energyDensity = out.energyDensity+e(beta)*r.energyDensity;
        out.energy = out.energy+e(beta)*r.energy;
        out.aliveEnergy = out.aliveEnergy+e(beta)*r.aliveEnergy;
    end
end
out.capturedEnergy = out.energy;
out.outsideObservationEnergy = max(0,sum(out.aliveEnergy,2)-sum(out.energy,2));
out.emittedEnergy = (sum(p)+sum(s))*double(first.t(:)>=bands.eventTime);
out.bandEnergy = [sum(p) sum(s)];
out.energyUnits = 'J';
out.energyDensityUnits = 'J per observation-cell measure';
out.unitInjectedEnergy = false;
out.composition = 'sum of physical weighted band energies times conditional unit responses';
end
