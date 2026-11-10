function [t, u] = freq2time(U, dt, eps)
%FREQ2TIME  Time series from spectra sampled at complex frequencies omega_j + 1i*eps.
%   [t, u] = lamb.wavefield.freq2time(U, dt, eps)
%
%   U   Nf-by-... spectra at omega_j = j*dw + 1i*eps, j = 0..Nf-1 (Nf = Nt/2 + 1,
%       dw = 2*pi/(Nt*dt)), i.e. the non-negative half up to Nyquist.
%   u   Nt-by-... real time series,
%         u(t) = exp(eps*t) * (1/2pi) int U(omega + 1i*eps) exp(-1i*omega*t) domega.
%
%   e^{-i omega t} convention: the synthesis kernel is exp(-1i*omega*t), which is MATLAB's
%   fft (not ifft). For a real signal, U(-omega + 1i*eps) = conj(U(omega + 1i*eps)).
sz = size(U);  Nf = sz(1);  Nt = 2*(Nf - 1);
dw = 2*pi/(Nt*dt);
U  = reshape(U, Nf, []);
U(Nf, :) = real(U(Nf, :));                                  % Nyquist bin must be real
Ufull = [U; conj(U(Nf-1:-1:2, :))];                         % negative frequencies
t = (0:Nt-1).'*dt;
u = real(fft(Ufull, [], 1))*dw/(2*pi).*exp(eps*t);
u = reshape(u, [Nt, sz(2:end)]);
end
