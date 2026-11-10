function K = interfaceLameKernels(ux, uz, dux, duz, k)
%INTERFACELAMEKERNELS  Closed-form Born kernels for delta-rho, -lambda, -mu.
%
%   K = interfaceLameKernels(ux, uz, dux, duz, k) returns function handles
%   K.rho_P(z), K.rho_Sv(z), K.rho_Sh(z), K.lam_P(z), K.mu_P(z), K.mu_Sv(z),
%   K.mu_Sh(z) for a mode travelling along +x with displacement
%
%       u = ( ux(z), 0, uz(z) ) exp(i k x)
%
%   in whatever frame the caller uses. ux, uz, dux, duz are handles.
%
%   These are Xu et al. (2021) Appendix B rewritten in the strains
%
%       exx = i k ux,   ezz = duz,   exz = (i k uz + dux)/2
%
%   which removes the (r1, r2) sign convention from the caller's hands. The
%   mapping is r1 = -exx/k, dr2/k = -ezz/k and dr1/k - r2 = -2i exz/k, so
%
%       lam_P  = |exx + ezz|^2 / k^2                       (B.1, B.15)
%       mu_P   = [3|exx|^2 + 3|ezz|^2 + 2Re(exx' ezz) + 4|exz|^2]/k^2   (B.12)
%       mu_Sv  = [|exx|^2/20 + 2|ezz|^2/15 - 2Re(exx' ezz)/15
%                 + 7|exz|^2/30]/k^2                        (B.13)
%       mu_Sh  = [|exx|^2/12 + |exz|^2/6]/k^2               (B.14, CORRECTED)
%
%   mu_Sh is the corrected kernel, not the printed B.14/B.22: the published
%   B.11 carries two wrong angular factors, which propagate into its average.
%   bornKernelsRayleigh.m holds the evidence; testBornKernels.m pins it.
%
%   Being closed form, these integrate over depth with quadrature instead of
%   a per-depth angular average, which is what makes them usable inside
%   scholteModel and stoneleyModel.
%
%   See also BORNKERNELS, BORNKERNELSRAYLEIGH.

exx = @(z) 1i*k*ux(z);
ezz = @(z) duz(z);
exz = @(z) 0.5*(1i*k*uz(z) + dux(z));

K.rho_P  = @(z) abs(ux(z)).^2 + abs(uz(z)).^2;
K.rho_Sv = @(z) abs(ux(z)).^2 + 4*abs(uz(z)).^2;
K.rho_Sh = @(z) abs(ux(z)).^2;

K.lam_P  = @(z) abs(exx(z) + ezz(z)).^2 / k^2;

K.mu_P   = @(z) (3*abs(exx(z)).^2 + 3*abs(ezz(z)).^2 ...
                 + 2*real(conj(exx(z)).*ezz(z)) + 4*abs(exz(z)).^2) / k^2;
K.mu_Sv  = @(z) (abs(exx(z)).^2/20 + 2*abs(ezz(z)).^2/15 ...
                 - 2*real(conj(exx(z)).*ezz(z))/15 + 7*abs(exz(z)).^2/30) / k^2;
K.mu_Sh  = @(z) (abs(exx(z)).^2/12 + abs(exz(z)).^2/6) / k^2;
end
