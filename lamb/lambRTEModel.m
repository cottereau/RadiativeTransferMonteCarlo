function M = lambRTEModel(varargin)
%LAMBRTEMODEL  Lamb (half-space Green's function) model of one planar interface for the MC.
%   M = lambRTEModel(iface, omega)
%   M = lambRTEModel(iface, omega, 'dsrc', d, 'zetaMin', 1e-3, 'NPhiBins', 360)
%
%   iface  lamb.Interface (lower solid at z > 0; upper solid, fluid or vacuum at z < 0)
%   omega  angular frequency of the MC run
%   dsrc   distance below which a source uses the Lamb emission K1 (default 3 lambda_S of
%          the lower solid). Beyond it: whole-space emission + Zoeppritz crossing.
%   zetaMin minimum k_beta1*|z| of MOMENT sources
%
%   Returns a LambModel object (see LAMBMODEL for its properties). Fields: .iface .omega .T (Lamb tables) .Tws (whole-space tables of the lower solid)
%   .ifaceFlip, .Tflip (solid-solid: the interface seen from the upper side, and its tables)
%   .dsrc .dslab .zetaMin
%   .kb1 (k_beta1) .phiBins .Rneg (9-by-9-by-NPhiBins, rotationMatrix(-phi))
%   .chan (struct array of channels: .name .side (1: z>0, 2: z<0) .mode (1 P, 2 S) .c)
%   .iw (struct array: .c speed) .leaky (.c speed, .rate per unit length, .sinf fluid angle)
%
%   See also INITIALIZEPARTICLESLAMB, PROPAGATEPARTICLELAMB, OBSERVELAMB.
M = LambModel(varargin{:});
end
