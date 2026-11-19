///|
/// Special functions the distributions need. Signatures only: each one
/// promises a domain and an accuracy bound in ULP, checked against an
/// oracle (SciPy, mpmath or Luna-Flow BinFloat) in a separate module.

///|
/// Natural log of |Γ(x)|.
pub fn lgamma(x : Double) -> Double {
  ignore(x)
  ...
}

///|
/// Regularized incomplete beta I_x(a, b), for a > 0, b > 0, 0 <= x <= 1.
/// The upper tail is a separate function so that callers never form
/// 1 - I_x(a, b) themselves.
pub fn betainc(a : Double, b : Double, x : Double) -> Double {
  ignore((a, b, x))
  ...
}

///|
/// 1 - I_x(a, b), computed without cancellation.
pub fn betaincc(a : Double, b : Double, x : Double) -> Double {
  ignore((a, b, x))
  ...
}

///|
/// Inverse of `betainc` in x.
pub fn betaincinv(a : Double, b : Double, p : Double) -> Double {
  ignore((a, b, p))
  ...
}

///|
/// Regularized lower incomplete gamma P(a, x).
pub fn gammainc(a : Double, x : Double) -> Double {
  ignore((a, x))
  ...
}

///|
/// Regularized upper incomplete gamma Q(a, x) = 1 - P(a, x).
pub fn gammaincc(a : Double, x : Double) -> Double {
  ignore((a, x))
  ...
}

///|
/// Standard normal lower tail Φ(x) and its log.
pub fn ndtr(x : Double) -> Double {
  ignore(x)
  ...
}

///|
pub fn log_ndtr(x : Double) -> Double {
  ignore(x)
  ...
}

///|
/// Φ⁻¹(p).
pub fn ndtri(p : Double) -> Double {
  ignore(p)
  ...
}
