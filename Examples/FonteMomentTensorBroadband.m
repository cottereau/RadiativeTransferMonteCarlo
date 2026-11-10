function example = FonteMomentTensorBroadband(mechanism,waveKind,opts)
%FONTEMOMENTTENSORBROADBAND Fonte DC-R/O/T, EX-R, CLVD or general tensor.
if nargin<1, mechanism='DC'; end
if nargin<2, waveKind='R'; end
if nargin<3, opts=struct; end
example=fonteBenchmarkCase(mechanism,waveKind,opts);
end
