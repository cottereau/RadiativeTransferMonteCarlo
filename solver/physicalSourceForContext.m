function src = physicalSourceForContext(source, d, acoustics, material)
% Rebuild from public fields so edited options never use a stale closure.
if strcmp(source.type,'acoustic')
    if d ~= 3 || ~acoustics
        error('Fonte:UnsupportedGeometry', 'Acoustic sources require 3D acoustics.');
    end
    for key = {'Q','c','rho'}
        if ~isfield(source,key{1}) || isempty(source.(key{1}))
            error('Fonte:MissingSourceField', 'Physical source requires %s.',key{1});
        end
    end
    if ~isfield(source,'F'), source.F = []; end
    src = sourceAcoustic(source.Q,source.F,source.c,source.rho,source);
    if nargin < 4, return; end
    checkMaterial(src.c,material.v,'c');
    if ~isempty(material.rho), checkMaterial(src.rho,material.rho,'rho'); end
    return
end
if d ~= 3 || acoustics
    error('Fonte:UnsupportedGeometry', 'Physical force/moment sources require 3D elasticity.');
end
required = {'vp','vs','rho'};
if strcmp(source.type,'point_force'), required{end+1} = 'F';
else, required{end+1} = 'M'; end
for i = 1:numel(required)
    if ~isfield(source,required{i}) || isempty(source.(required{i}))
        error('Fonte:MissingSourceField', 'Physical source requires %s.',required{i});
    end
end
if strcmp(source.type,'point_force')
    src = sourcePointForce(source.F,source.vp,source.vs,source.rho,source);
else
    src = sourceMomentTensor(source.M,source.vp,source.vs,source.rho,source);
end
if nargin < 4, return; end
for name = {'vp','vs','rho'}
    key = name{1};
    value = material.(key);
    % Legacy materials may omit absolute density. It affects later physical
    % energy, not the normalized trajectories; check it when available.
    if strcmp(key,'rho') && isempty(value), continue; end
    checkMaterial(src.(key),value,key);
end
end

function checkMaterial(sourceValue,value,key)
validateattributes(value,{'numeric'},{'real','finite','scalar','positive'});
if abs(sourceValue-value) > 1e-10*max(abs(sourceValue),abs(value))
    error('Fonte:MaterialMismatch', 'Source %s differs from the transport material.',key);
end
end
