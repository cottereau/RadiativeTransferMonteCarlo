function h = heaviside(x)

  h = zeros(size(x));
  h(x > 0) = 1.0;
  h(x == 0) = 0.;