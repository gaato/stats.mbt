///|
/// The rest of the special functions the distributions and tests need.
/// Error policy, after Boost.Math's: these plain functions return NaN
/// outside their domain (documented per function) so that hot loops stay
/// free of error handling; the distribution constructors check parameters,
/// so a NaN here means a caller bypassed them.

///|
pub fn erf(x : Double) -> Double {
  ignore(x)
  ...
}

///|
/// 1 - erf(x) without cancellation for large x.
pub fn erfc(x : Double) -> Double {
  ignore(x)
  ...
}

///|
pub fn erfinv(p : Double) -> Double {
  ignore(p)
  ...
}

///|
/// ψ(x) = Γ'(x)/Γ(x); fitting the gamma shape solves ψ(a) - ln a = c.
pub fn digamma(x : Double) -> Double {
  ignore(x)
  ...
}

///|
pub fn trigamma(x : Double) -> Double {
  ignore(x)
  ...
}

///|
/// ln B(a, b) = ln Γ(a) + ln Γ(b) - ln Γ(a + b), without the cancellation
/// of the three-term form when one argument is large.
pub fn lbeta(a : Double, b : Double) -> Double {
  ignore((a, b))
  ...
}

///|
/// ln C(n, k), for binomial and hypergeometric masses.
pub fn lchoose(n : Double, k : Double) -> Double {
  ignore((n, k))
  ...
}

///|
/// Inverse of the upper incomplete gamma in x.
pub fn gammainccinv(a : Double, q : Double) -> Double {
  ignore((a, q))
  ...
}

///|
/// Inverse of the lower incomplete gamma in x.
pub fn gammaincinv(a : Double, p : Double) -> Double {
  ignore((a, p))
  ...
}

///|
/// Owen's T function, for the bivariate normal and skew-normal CDFs.
pub fn owens_t(h : Double, a : Double) -> Double {
  ignore((h, a))
  ...
}

///|
/// Kolmogorov's limiting distribution, for the KS test p-value.
pub fn kolmogorov_sf(x : Double) -> Double {
  ignore(x)
  ...
}

///|
/// Studentized range distribution, for Tukey's HSD.
pub fn studentized_range_sf(q : Double, k : Double, df : Double) -> Double {
  ignore((q, k, df))
  ...
}
