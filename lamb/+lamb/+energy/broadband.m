function EPb = broadband(iface, src, omegaVec, varargin)
%BROADBAND  Wavelet-weighted partition: W_m = int |s(omega)|^2 W_m(omega) domega.
%   EPb = lamb.energy.broadband(iface, src, omegaVec, <energyPartition options>)
%   EPb.perOmega : struct array of lamb.energyPartition outputs
%   EPb.W        : integrated powers (.body struct, .interface, .leaky, .total)
%   EPb.pct      : percentages of EPb.W.total
omegaVec = omegaVec(:).';
for i = numel(omegaVec):-1:1
    EPs(i) = lamb.energyPartition(iface, src, omegaVec(i), varargin{:});
end
wt = abs(src.wavelet(omegaVec)).^2;
I  = @(v) trapz(omegaVec, wt.*v);
for n = fieldnames(EPs(1).W.body).'
    EPb.W.body.(n{1}) = I(arrayfun(@(e) e.W.body.(n{1}), EPs));
end
EPb.W.interface = I(arrayfun(@(e) e.W.ifTotal, EPs));
EPb.W.leaky     = I(arrayfun(@(e) sum([e.W.leaky.W]), EPs));
EPb.W.total     = I(arrayfun(@(e) e.W.total, EPs));
EPb.pct.body      = structfun(@(v) 100*v/EPb.W.total, EPb.W.body, 'UniformOutput', false);
EPb.pct.interface = 100*EPb.W.interface/EPb.W.total;
EPb.pct.leaky     = 100*EPb.W.leaky/EPb.W.total;
EPb.perOmega = EPs;
end
