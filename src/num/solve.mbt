///|
/// Iterative numerics with an honest status. Every solver and series
/// evaluator returns how it ended, so a caller (a quantile function, a
/// special function) can refuse to report a value it could not justify.

///|
#derive(Eq, Debug)
pub(all) enum Status {
  Converged
  MaxIterations
  /// The function did not change sign over the bracket.
  NoBracket
  /// A NaN appeared during the iteration.
  NaN
}

///|
pub struct Solution {
  value : Double
  status : Status
  iterations : Int
  /// Width of the final bracket, or the last step for open methods.
  error_bound : Double
}

///|
/// Root of a function known to change sign on [lo, hi]. Bracketing
/// methods cannot diverge, which is what quantile functions need: every
/// distribution's quantile without a closed form is a root of
/// F(x) - p on a bracket. Uses Brent's method (inverse quadratic
/// interpolation with bisection as the fallback).
pub fn find_root(
  f : (Double) -> Double,
  lo : Double,
  hi : Double,
  x_tol? : Double = 0.0,
  max_iterations? : Int = 200,
) -> Solution {
  ignore((f, lo, hi, x_tol, max_iterations))
  ...
}

///|
/// Root of a monotone function with its derivative, Newton steps kept
/// inside a bracket (falling back to bisection when a step leaves it).
/// For quantiles when both the distribution function and the density are
/// available: quadratic convergence without the risk of open Newton.
pub fn find_root_newton(
  f : (Double) -> (Double, Double),
  lo : Double,
  hi : Double,
  start : Double,
  x_tol? : Double = 0.0,
  max_iterations? : Int = 100,
) -> Solution {
  ignore((f, lo, hi, start, x_tol, max_iterations))
  ...
}

///|
/// Continued fraction b0 + a1/(b1 + a2/(b2 + ...)) by the modified Lentz
/// method; `term(n)` returns (aₙ, bₙ) for n >= 1. The incomplete beta and
/// gamma functions are evaluated this way in the region where their power
/// series converge slowly.
pub fn continued_fraction(
  b0 : Double,
  term : (Int) -> (Double, Double),
  rel_tol? : Double = 1.0e-16,
  max_terms? : Int = 1000,
) -> Solution {
  ignore((b0, term, rel_tol, max_terms))
  ...
}

///|
/// Σ term(n) for n >= 0 until a term falls below rel_tol of the sum,
/// compensated.
pub fn sum_series(
  term : (Int) -> Double,
  rel_tol? : Double = 1.0e-16,
  max_terms? : Int = 10000,
) -> Solution {
  ignore((term, rel_tol, max_terms))
  ...
}

///|
/// ∫ f over [a, b], adaptive Gauss–Kronrod, for moments and probabilities
/// of distributions with a density but no closed-form distribution
/// function.
pub fn integrate(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  abs_tol? : Double = 0.0,
  rel_tol? : Double = 1.0e-12,
) -> Solution {
  ignore((f, a, b, abs_tol, rel_tol))
  ...
}

///|
/// Neumaier's compensated sum: robust also when a later term is larger
/// than the running sum. A value type, so it can be folded or merged.
pub struct CompensatedSum {
  sum : Double
  c : Double
}

///|
pub fn CompensatedSum::zero() -> CompensatedSum {
  { sum: 0.0, c: 0.0 }
}

///|
pub fn CompensatedSum::add(self : CompensatedSum, x : Double) -> CompensatedSum {
  let t = self.sum + x
  let c = if self.sum.abs() >= x.abs() {
    self.c + (self.sum - t + x)
  } else {
    self.c + (x - t + self.sum)
  }
  { sum: t, c }
}

///|
pub fn CompensatedSum::merge(
  self : CompensatedSum,
  other : CompensatedSum,
) -> CompensatedSum {
  let s = self.add(other.sum)
  { sum: s.sum, c: s.c + other.c }
}

///|
pub fn CompensatedSum::value(self : CompensatedSum) -> Double {
  self.sum + self.c
}
