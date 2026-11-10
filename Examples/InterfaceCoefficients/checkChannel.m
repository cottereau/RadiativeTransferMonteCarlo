function ok = checkChannel(label, theory, mc, tol)
%CHECKCHANNEL  Print and pass/fail one Monte Carlo channel against theory.
%
%   ok = checkChannel(label, theory, mc, tol) prints theory vs the
%   Monte Carlo fraction (mc) for a named channel, flags PASS/FAIL against
%   an absolute tolerance tol (loose enough to absorb sampling noise at
%   the particle counts used in this folder's examples, not floating-point
%   noise), and returns whether it passed.
%
%   Shared by every script in Examples/InterfaceCoefficients/.

ok = abs(theory-mc) < tol;
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('    %-10s theory %7.4f   MC %7.4f   |diff| %7.4f   [%s]\n', ...
    label, theory, mc, abs(theory-mc), tag);
end
