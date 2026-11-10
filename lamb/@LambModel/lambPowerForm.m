function Qt = lambPowerForm(M, z)
%LAMBPOWERFORM  Total radiated-power form at depth z (physical frame): W(s) = real(s'*Qt*s).
%   Same tables and mirror as lambScatteringModel's powerAt (upper solid: M.Tflip, mirrored).
w = M.omega;
if z < 0 && strcmp(M.iface.kind, 'solid-solid')
    Q = lamb.rte.channelQ(M.Tflip, -(w/M.iface.upper.beta)*z, w);  D = diag([1 1 -1 1 1 1 1 -1 -1]);
else
    Q = lamb.rte.channelQ(M.T, M.kb1*z, w);  D = eye(9);
end
Qt = zeros(9);
for f = fieldnames(Q.body).', Qt = Qt + Q.body.(f{1}); end
for j = 1:numel(Q.iw), Qt = Qt + Q.iw{j}; end
Qt = D*Qt*D;
end
