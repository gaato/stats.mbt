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
  iterations : Int32
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
  x_tol~ : Double = 0.0,
  max_iterations~ : Int32 = 200,
) -> Solution {
  let mut a = lo
  let mut b = hi
  let mut fa = f(a)
  let mut fb = f(b)
  if fa.is_nan() || fb.is_nan() {
    return {
      value: 0.0 / 0.0,
      status: NaN,
      iterations: 0,
      error_bound: 1.0 / 0.0,
    }
  }
  if fa == 0.0 {
    return { value: a, status: Converged, iterations: 0, error_bound: 0.0, }
  }
  if fb == 0.0 {
    return { value: b, status: Converged, iterations: 0, error_bound: 0.0, }
  }
  if (fa > 0.0) == (fb > 0.0) {
    return {
      value: 0.0 / 0.0,
      status: NoBracket,
      iterations: 0,
      error_bound: 1.0 / 0.0,
    }
  }
  // Brent (1973): keep b the best estimate, c the previous b, and a bracket
  // [b, c] with a sign change; interpolate when it shrinks fast enough.
  let mut c = a
  let mut fc = fa
  let mut d = b - a
  let mut e = d
  for it in 1..<=max_iterations {
    if (fb > 0.0) == (fc > 0.0) {
      c = a
      fc = fa
      d = b - a
      e = d
    }
    if fc.abs() < fb.abs() {
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
    }
    let tol = 2.0 * 2.220446049250313e-16 * b.abs() + 0.5 * x_tol
    let m = 0.5 * (c - b)
    if m.abs() <= tol || fb == 0.0 {
      return {
        value: b,
        status: Converged,
        iterations: it,
        error_bound: m.abs() * 2.0,
      }
    }
    if e.abs() >= tol && fa.abs() > fb.abs() {
      let s0 = fb / fa
      let mut pp = 0.0
      let mut qq = 0.0
      if a == c {
        pp = 2.0 * m * s0
        qq = 1.0 - s0
      } else {
        let q0 = fa / fc
        let r = fb / fc
        pp = s0 * (2.0 * m * q0 * (q0 - r) - (b - a) * (r - 1.0))
        qq = (q0 - 1.0) * (r - 1.0) * (s0 - 1.0)
      }
      if pp > 0.0 {
        qq = -qq
      } else {
        pp = -pp
      }
      if 2.0 * pp < 3.0 * m * qq - (tol * qq).abs() && pp < (0.5 * e * qq).abs() {
        e = d
        d = pp / qq
      } else {
        d = m
        e = m
      }
    } else {
      d = m
      e = m
    }
    a = b
    fa = fb
    b = b + (if d.abs() > tol { d } else if m > 0.0 { tol } else { -tol })
    fb = f(b)
    if fb.is_nan() {
      return { value: b, status: NaN, iterations: it, error_bound: 1.0 / 0.0, }
    }
  }
  {
    value: b,
    status: MaxIterations,
    iterations: max_iterations,
    error_bound: (c - b).abs(),
  }
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
  x_tol~ : Double = 0.0,
  max_iterations~ : Int32 = 100,
) -> Solution {
  let mut lo = lo
  let mut hi = hi
  let (flo, _) = f(lo)
  let (fhi, _) = f(hi)
  if flo.is_nan() || fhi.is_nan() {
    return {
      value: 0.0 / 0.0,
      status: NaN,
      iterations: 0,
      error_bound: 1.0 / 0.0,
    }
  }
  if (flo > 0.0) == (fhi > 0.0) && flo != 0.0 && fhi != 0.0 {
    return {
      value: 0.0 / 0.0,
      status: NoBracket,
      iterations: 0,
      error_bound: 1.0 / 0.0,
    }
  }
  let increasing = flo < fhi
  let mut x = if start > lo && start < hi { start } else { 0.5 * (lo + hi) }
  for it in 1..<=max_iterations {
    let (fx, dfx) = f(x)
    if fx.is_nan() {
      return { value: x, status: NaN, iterations: it, error_bound: 1.0 / 0.0, }
    }
    if fx == 0.0 {
      return { value: x, status: Converged, iterations: it, error_bound: 0.0, }
    }
    // shrink the bracket with the sign of f
    if (fx > 0.0) == increasing {
      hi = x
    } else {
      lo = x
    }
    let newton = if dfx != 0.0 && !dfx.is_nan() {
      x - fx / dfx
    } else {
      0.0 / 0.0
    }
    let mid = if lo > 0.0 && hi > 4.0 * lo {
      // the bracket spans decades: halve it in the logarithm
      (lo * hi).sqrt()
    } else {
      0.5 * (lo + hi)
    }
    let next = if newton > lo && newton < hi { newton } else { mid }
    let step = (next - x).abs()
    x = next
    if step <= x_tol + 4.440892098500626e-16 * x.abs() ||
      hi - lo <= x_tol + 4.440892098500626e-16 * x.abs() {
      return { value: x, status: Converged, iterations: it, error_bound: step, }
    }
  }
  {
    value: x,
    status: MaxIterations,
    iterations: max_iterations,
    error_bound: hi - lo,
  }
}

///|
/// Continued fraction b0 + a1/(b1 + a2/(b2 + ...)) by the modified Lentz
/// method; `term(n)` returns (aₙ, bₙ) for n >= 1. The incomplete beta and
/// gamma functions are evaluated this way in the region where their power
/// series converge slowly.
pub fn continued_fraction(
  b0 : Double,
  term : (Int32) -> (Double, Double),
  rel_tol~ : Double = 1.0e-16,
  max_terms~ : Int32 = 1000,
) -> Solution {
  // Modified Lentz (Thompson and Barnett 1986).
  let tiny = 1.0e-300
  let mut f = if b0 == 0.0 { tiny } else { b0 }
  let mut c = f
  let mut d = 0.0
  for n in 1..<=max_terms {
    let (an, bn) = term(n)
    d = bn + an * d
    if d == 0.0 {
      d = tiny
    }
    c = bn + an / c
    if c == 0.0 {
      c = tiny
    }
    d = 1.0 / d
    let delta = c * d
    f = f * delta
    if f.is_nan() {
      return { value: f, status: NaN, iterations: n, error_bound: 1.0 / 0.0, }
    }
    if (delta - 1.0).abs() <= rel_tol {
      return {
        value: f,
        status: Converged,
        iterations: n,
        error_bound: f.abs() * rel_tol,
      }
    }
  }
  {
    value: f,
    status: MaxIterations,
    iterations: max_terms,
    error_bound: 1.0 / 0.0,
  }
}

///|
/// Σ term(n) for n >= 0 until a term falls below rel_tol of the sum,
/// compensated.
pub fn sum_series(
  term : (Int32) -> Double,
  rel_tol~ : Double = 1.0e-16,
  max_terms~ : Int32 = 10000,
) -> Solution {
  let mut acc = CompensatedSum::zero()
  for n in 0..<max_terms {
    let t = term(n)
    if t.is_nan() {
      return { value: t, status: NaN, iterations: n, error_bound: 1.0 / 0.0, }
    }
    acc = acc.add(t)
    if t.abs() <= rel_tol * acc.value().abs() {
      return {
        value: acc.value(),
        status: Converged,
        iterations: n + 1,
        error_bound: t.abs(),
      }
    }
  }
  {
    value: acc.value(),
    status: MaxIterations,
    iterations: max_terms,
    error_bound: 1.0 / 0.0,
  }
}

///|
/// ∫ f over [a, b], adaptive Gauss–Kronrod, for moments and probabilities
/// of distributions with a density but no closed-form distribution
/// function.
pub fn integrate(
  f : (Double) -> Double,
  a : Double,
  b : Double,
  abs_tol~ : Double = 0.0,
  rel_tol~ : Double = 1.0e-12,
) -> Solution {
  // Double exponential (tanh-sinh) quadrature of Takahasi and Mori (1974):
  // x = c + h tanh(π/2 sinh t); the error falls exponentially in the number
  // of points for analytic integrands, and endpoint singularities are
  // harmless because no node sits on an endpoint.
  let c = 0.5 * (a + b)
  let h = 0.5 * (b - a)
  let half_pi = 0.5 * @math.PI
  let node = (t : Double) => {
    let u = half_pi * @math.sinh(t)
    let ch = @math.cosh(u)
    let w = half_pi * @math.cosh(t) / (ch * ch)
    // distance from the nearer endpoint without cancellation: 1 - tanh|u|
    let gap = 2.0 / (1.0 + @math.exp(2.0 * u.abs()))
    let x = if u >= 0.0 { b - h * gap } else { a + h * gap }
    (x, w)
  }
  let t_max = 3.2
  let mut step = 1.0
  let mut total = f(c) * half_pi
  let mut k = 1
  while k.to_double() * step <= t_max {
    let (x1, w1) = node(k.to_double() * step)
    let (x2, w2) = node(-k.to_double() * step)
    total = total + w1 * f(x1) + w2 * f(x2)
    k = k + 1
  }
  let mut estimate = h * step * total
  for level in 1..<=12 {
    step = step * 0.5
    // only the new (odd) nodes are evaluated at each level
    let mut k = 1
    while k.to_double() * step <= t_max {
      let (x1, w1) = node(k.to_double() * step)
      let (x2, w2) = node(-k.to_double() * step)
      total = total + w1 * f(x1) + w2 * f(x2)
      k = k + 2
    }
    let next = h * step * total
    let err = (next - estimate).abs()
    estimate = next
    if estimate.is_nan() {
      return {
        value: estimate,
        status: NaN,
        iterations: level,
        error_bound: 1.0 / 0.0,
      }
    }
    if level >= 3 && err <= abs_tol + rel_tol * estimate.abs() {
      return {
        value: estimate,
        status: Converged,
        iterations: level,
        error_bound: err,
      }
    }
  }
  {
    value: estimate,
    status: MaxIterations,
    iterations: 12,
    error_bound: 1.0 / 0.0,
  }
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
  { sum: 0.0, c: 0.0, }
}

///|
pub fn CompensatedSum::add(self : CompensatedSum, x : Double) -> CompensatedSum {
  let t = self.sum + x
  let c = if self.sum.abs() >= x.abs() {
    self.c + (self.sum - t + x)
  } else {
    self.c + (x - t + self.sum)
  }
  { sum: t, c, }
}

///|
pub fn CompensatedSum::merge(
  self : CompensatedSum,
  other : CompensatedSum,
) -> CompensatedSum {
  let s = self.add(other.sum)
  { sum: s.sum, c: s.c + other.c, }
}

///|
pub fn CompensatedSum::value(self : CompensatedSum) -> Double {
  self.sum + self.c
}

///|
/// A domain error for the checked variants of special functions; the
/// plain variants return NaN instead.
#non_exhaustive
#derive(Debug)
pub(all) suberror DomainError {
  OutsideDomain(function~ : String, argument~ : Double)
  Overflow(function~ : String)
}
