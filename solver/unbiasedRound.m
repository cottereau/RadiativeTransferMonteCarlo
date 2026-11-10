function k = unbiasedRound(x)
% Randomized rounding with E[k] = x, so small mixed P/S launches carry the
% right modal energy on average (round(1/17*4) = 0 would drop all P).
% Values within 1e-9 of an integer are taken exactly (no RNG draw).
k = round(x);
if abs(x-k) <= 1e-9*max(1,abs(x)), return; end
k = floor(x); k = k+(rand < x-k);
end
