function example = FontePointForceBroadband(waveKind,opts)
%FONTEPOINTFORCEBROADBAND Fonte F-R, F-O or F-T configuration and run.
if nargin<1, waveKind='R'; end
if nargin<2, opts=struct; end
example=fonteBenchmarkCase('F',waveKind,opts);
end
