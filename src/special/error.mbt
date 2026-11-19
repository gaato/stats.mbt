///|
/// The error function family and the normal distribution's tails.
///
/// erf uses its all-positive series below 1 (no cancellation); erfc above 1
/// is e^{-x²} times the scaled function erfcx, a continued fraction that
/// never forms e^{-x²} itself, and e^{-x²} is computed from an exact split
/// of x² so that the large exponent does not cost digits.

///|
let inv_sqrt_pi : Double = 0.5641895835477563

///|
let sqrt_2 : Double = 1.4142135623730951

///|
let ln_sqrt_2pi : Double = 0.9189385332046728

///|
/// e^{-x²/2·scale·2} without rounding x²: x = hi + lo with hi on a 1/64 grid,
/// so hi² is exact and lo(x + hi) is small. `half` gives e^{-x²/2}.
fn exp_neg_sq(x : Double, half : Bool) -> Double {
  let ax = x.abs()
  if ax > 1.0e6 {
    return 0.0
  }
  let hi = @math.round(ax * 64.0) / 64.0
  let lo = ax - hi
  let k = if half { 0.5 } else { 1.0 }
  @math.exp(-k * hi * hi) * @math.exp(-k * lo * (ax + hi))
}

///|
/// erf(x) = (2x/√π) e^{-x²} Σ (2x²)^n / (1·3·...·(2n+1)); every term is
/// positive.
fn erf_series(x : Double) -> Double {
  let x2 = 2.0 * x * x
  let mut term = 1.0
  let mut sum = 1.0
  let mut n = 0
  while n < 200 {
    n = n + 1
    term = term * x2 / (2 * n + 1).to_double()
    sum = sum + term
    if term <= 1.0e-17 * sum {
      break
    }
  }
  2.0 * x * inv_sqrt_pi * exp_neg_sq(x, false) * sum
}

///|
/// erfcx(x) = e^{x²} erfc(x) for x >= 1 by the continued fraction
/// 1/(√π (x + (1/2)/(x + 1/(x + (3/2)/(x + ...))))).
fn erfcx_cf(x : Double) -> Double {
  if x > 1.0e8 {
    return inv_sqrt_pi / x
  }
  let r = @num.continued_fraction(
    x,
    n => (0.5 * n.to_double(), x),
    max_terms=2000,
  )
  inv_sqrt_pi / r.value
}

///|
/// The error function. Relative error near 1e-16 everywhere.
pub fn erf(x : Double) -> Double {
  if x.is_nan() {
    return x
  }
  let ax = x.abs()
  let v = if ax < 1.0 {
    erf_series(ax)
  } else if ax > 6.0 {
    1.0
  } else {
    1.0 - exp_neg_sq(ax, false) * erfcx_cf(ax)
  }
  if x < 0.0 {
    -v
  } else {
    v
  }
}

///|
/// 1 - erf(x), with full relative accuracy in the upper tail (down to the
/// underflow near x = 27.2).
pub fn erfc(x : Double) -> Double {
  if x.is_nan() {
    return x
  }
  if x < 1.0 {
    if x > -1.0 {
      1.0 - erf(x)
    } else {
      2.0 - erfc(-x)
    }
  } else {
    exp_neg_sq(x, false) * erfcx_cf(x)
  }
}

///|
/// The scaled complementary error function e^{x²} erfc(x), finite for all
/// x >= 0 where erfc itself underflows.
pub fn erfcx(x : Double) -> Double {
  if x.is_nan() {
    return x
  }
  if x >= 1.0 {
    erfcx_cf(x)
  } else {
    // e^{x²} is moderate here; erfc(x) has no cancellation for x < 1
    erfc(x) * @math.exp(x * x)
  }
}

///|
/// Standard normal lower tail Φ(z). In the tails the e^{-z²/2} factor is
/// formed from z itself, not from z/√2, so no rounding of the argument is
/// amplified by the exponent.
pub fn ndtr(z : Double) -> Double {
  if z.is_nan() {
    return z
  }
  if z < -sqrt_2 {
    0.5 * exp_neg_sq(z, true) * erfcx_cf(-z / sqrt_2)
  } else if z > sqrt_2 {
    1.0 - 0.5 * exp_neg_sq(z, true) * erfcx_cf(z / sqrt_2)
  } else {
    0.5 * erfc(-z / sqrt_2)
  }
}

///|
/// ln Φ(z), finite far below where Φ underflows.
pub fn log_ndtr(z : Double) -> Double {
  if z.is_nan() {
    return z
  }
  if z < -sqrt_2 {
    -0.5 * z * z + @math.ln(0.5 * erfcx_cf(-z / sqrt_2))
  } else if z > 0.0 {
    // ln Φ(z) is near 0 here; ln1p of the small lower tail keeps its digits
    @math.ln_1p(-ndtr(-z))
  } else {
    @math.ln(ndtr(z))
  }
}

///|
/// The standard normal density.
fn npdf(z : Double) -> Double {
  exp_neg_sq(z, true) * 0.3989422804014327
}

///|
/// Φ⁻¹(p). A rational start (Abramowitz and Stegun 26.2.23, error below
/// 4.5e-4) refined by Newton steps on ln Φ, which keeps full relative
/// accuracy down to the smallest positive p.
pub fn ndtri(p : Double) -> Double {
  if p.is_nan() || p < 0.0 || p > 1.0 {
    return 0.0 / 0.0
  }
  if p == 0.0 {
    return -1.0 / 0.0
  }
  if p == 1.0 {
    return 1.0 / 0.0
  }
  if p == 0.5 {
    return 0.0
  }
  // solve in the lower tail; 1 - p is exact for p >= 1/2
  let (q, sign) = if p < 0.5 { (p, -1.0) } else { (1.0 - p, 1.0) }
  let t = (-2.0 * @math.ln(q)).sqrt()
  let mut x = -(t -
    (2.515517 + 0.802853 * t + 0.010328 * t * t) /
    (1.0 + 1.432788 * t + 0.189269 * t * t + 0.001308 * t * t * t))
  let lq = @math.ln(q)
  for _ in 0..<10 {
    let f = log_ndtr(x) - lq
    // d/dx ln Φ(x) = φ(x)/Φ(x)
    let phi = ndtr(x)
    let g = npdf(x) / phi
    let step = f / g
    x = x - step
    if step.abs() <= 1.0e-16 * x.abs() {
      break
    }
  }
  sign * -x
}

///|
/// erf⁻¹(y) for -1 < y < 1: a start from Φ⁻¹, then Newton steps on erf
/// near 0 and on erfc near ±1, so the result keeps its relative accuracy at
/// both ends.
pub fn erfinv(y : Double) -> Double {
  if y.is_nan() || y < -1.0 || y > 1.0 {
    return 0.0 / 0.0
  }
  if y == 1.0 {
    return 1.0 / 0.0
  }
  if y == -1.0 {
    return -1.0 / 0.0
  }
  if y == 0.0 {
    return y
  }
  let ay = y.abs()
  let c = 1.0 - ay
  let mut x = -ndtri(0.5 * c) / sqrt_2
  for _ in 0..<10 {
    let d = 2.0 * inv_sqrt_pi * exp_neg_sq(x, false)
    let step = if ay < 0.5 { (erf(x) - ay) / d } else { (c - erfc(x)) / d }
    x = x - step
    if step.abs() <= 1.0e-16 * x.abs() {
      break
    }
  }
  if y < 0.0 {
    -x
  } else {
    x
  }
}
