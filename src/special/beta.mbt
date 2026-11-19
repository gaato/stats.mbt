///|
/// The regularised incomplete beta function, its inverse, and the
/// remaining special functions (Owen's T, the Kolmogorov distribution).
///
/// x^a (1-x)^b / B(a, b) is formed, when both a and b are large, as
/// √(ab/(2π(a+b))) exp(a·ln1pmx(d/a) + b·ln1pmx(-d/b) + δ(a+b) - δ(a) - δ(b))
/// with d = xb - (1-x)a, which keeps its digits where the direct logarithm
/// would subtract numbers of size a ln a.

///|
/// x^a y^b / B(a, b) with y = 1 - x given separately.
fn beta_prefix(a : Double, b : Double, x : Double, y : Double) -> Double {
  if a >= 8.0 && b >= 8.0 {
    let d = x * b - y * a
    let lx = if x < 0.5 { @math.ln(x) } else { @math.ln_1p(-y) }
    let ly = if y < 0.5 { @math.ln(y) } else { @math.ln_1p(-x) }
    // ln(1+u) - u with u = d/a; far from the mean 1 + u = x(a+b)/a is tiny
    // and ln1p(u) would cancel, so take ln x + ln(1 + b/a) instead
    let u = d / a
    let v = -d / b
    let ta = if u.abs() < 0.5 {
      log1pmx(u)
    } else {
      lx + @math.ln_1p(b / a) - u
    }
    let tb = if v.abs() < 0.5 {
      log1pmx(v)
    } else {
      ly + @math.ln_1p(a / b) - v
    }
    (a * b / (2.0 * @math.PI * (a + b))).sqrt() *
    @math.exp(
      a * ta +
      b * tb +
      stirling_delta(a + b) -
      stirling_delta(a) -
      stirling_delta(b),
    )
  } else {
    let lx = if x < 0.5 { @math.ln(x) } else { @math.ln_1p(-y) }
    let ly = if y < 0.5 { @math.ln(y) } else { @math.ln_1p(-x) }
    @math.exp(a * lx + b * ly - lbeta(a, b))
  }
}

///|
/// I_x(a, b) by the continued fraction (DLMF 8.17.22), valid and fast for
/// x < (a + 1)/(a + b + 2).
fn ibeta_cf(a : Double, b : Double, x : Double, y : Double) -> Double {
  let r = @num.continued_fraction(
    1.0,
    n => {
      let m = (n / 2).to_double()
      let d = if n % 2 == 1 {
        -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
      } else {
        m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
      }
      (d, 1.0)
    },
    max_terms=100000,
  )
  beta_prefix(a, b, x, y) / (a * r.value)
}

///|
/// I_x(a, b) by its power series x^a/B(a,b) Σ (1-b)_n x^n / (n! (a+n)).
/// Every term is positive for b < 1, so it gives the small tail directly
/// where the continued fraction's side is near 1 and its complement would
/// cancel. `None` if it does not converge within the term limit.
fn ibeta_series(a : Double, b : Double, x : Double, y : Double) -> Double? {
  let ly = if y < 0.5 { @math.ln(y) } else { @math.ln_1p(-x) }
  let pre = beta_prefix(a, b, x, y) / @math.exp(b * ly)
  let mut c = 1.0
  let mut sum = 1.0 / a
  for n in 1..<=100000 {
    let nd = n.to_double()
    c = c * (nd - b) * x / nd
    let t = c / (a + nd)
    sum = sum + t
    if t <= 1.0e-17 * sum {
      return Some(pre * sum)
    }
  }
  None
}

///|
/// (I_x(a, b), 1 - I_x(a, b)) with y = 1 - x given separately, so that a
/// caller who knows 1 - x exactly (a t distribution's x²/(ν + x²)) keeps it.
/// The continued fraction is evaluated on whichever side converges, and the
/// complement taken on the other.
pub fn betainc_pair(
  a : Double,
  b : Double,
  x : Double,
  y : Double,
) -> (Double, Double) {
  if a.is_nan() || b.is_nan() || x.is_nan() || a <= 0.0 || b <= 0.0 {
    return (0.0 / 0.0, 0.0 / 0.0)
  }
  if x <= 0.0 {
    return (0.0, 1.0)
  }
  if y <= 0.0 {
    return (1.0, 0.0)
  }
  if x < (a + 1.0) / (a + b + 2.0) {
    let w = ibeta_cf(a, b, x, y)
    // a small a piles the mass near 0: w is near 1 and 1 - w cancels
    if w > 0.5 && a < 1.0 {
      match ibeta_series(b, a, y, x) {
        Some(c) => return (w, c)
        None => ()
      }
    }
    (w, 1.0 - w)
  } else {
    let w = ibeta_cf(b, a, y, x)
    if w > 0.5 && b < 1.0 {
      match ibeta_series(a, b, x, y) {
        Some(i) => return (i, w)
        None => ()
      }
    }
    (1.0 - w, w)
  }
}

///|
/// The regularised incomplete beta I_x(a, b), a, b > 0, 0 <= x <= 1.
pub fn betainc(a : Double, b : Double, x : Double) -> Double {
  betainc_pair(a, b, x, 1.0 - x).0
}

///|
/// 1 - I_x(a, b), computed without cancellation.
pub fn betaincc(a : Double, b : Double, x : Double) -> Double {
  betainc_pair(a, b, x, 1.0 - x).1
}

///|
/// x with I_x(a, b) = p: a start from Abramowitz and Stegun 26.5.22 (or
/// the power-law tails when a or b is below 1), then Newton steps on the
/// logarithm of the smaller tail inside [0, 1].
pub fn betaincinv(a : Double, b : Double, p : Double) -> Double {
  if a.is_nan() ||
    b.is_nan() ||
    p.is_nan() ||
    a <= 0.0 ||
    b <= 0.0 ||
    p < 0.0 ||
    p > 1.0 {
    return 0.0 / 0.0
  }
  if p == 0.0 {
    return 0.0
  }
  if p == 1.0 {
    return 1.0
  }
  // solve I_x(a, b) = p if p <= 1/2, else 1 - I_x(a, b) = 1 - p
  let upper = p > 0.5
  let target = if upper { 1.0 - p } else { p }
  let mut x0 = if a >= 1.0 && b >= 1.0 {
    let pp = if p < 0.5 { p } else { 1.0 - p }
    let t = (-2.0 * @math.ln(pp)).sqrt()
    let mut z = (2.30753 + 0.27061 * t) / (1.0 + t * (0.99229 + 0.04481 * t)) -
      t
    if p < 0.5 {
      z = -z
    }
    let al = (z * z - 3.0) / 6.0
    let h = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0))
    let w = z * (al + h).sqrt() / h -
      (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0)) *
      (al + 5.0 / 6.0 - 2.0 / (3.0 * h))
    a / (a + b * @math.exp(2.0 * w))
  } else {
    let lna = @math.ln(a / (a + b))
    let lnb = @math.ln(b / (a + b))
    let t = @math.exp(a * lna) / a
    let u = @math.exp(b * lnb) / b
    let w = t + u
    if p < t / w {
      @math.pow(a * w * p, 1.0 / a)
    } else {
      1.0 - @math.pow(b * w * (1.0 - p), 1.0 / b)
    }
  }
  if !(x0 > 0.0 && x0 < 1.0) {
    x0 = 0.5
  }
  // Work in the variable that is small at the solution: x for the lower
  // tail, y = 1 - x for the upper, and in u = ln of it for the Newton steps.
  let lt = @math.ln(target)
  let at = (s : Double) => if upper { (1.0 - s, s) } else { (s, 1.0 - s) }
  let value = (s : Double) => {
    let (x, y) = at(s)
    let (i, ic) = betainc_pair(a, b, x, y)
    if upper {
      ic
    } else {
      i
    }
  }
  // d/ds of the tail in s, which is positive for the lower tail in x and
  // negative for the upper tail in y... both tails decrease as their own
  // variable shrinks, so the derivative in s is the density for both
  let density = (s : Double) => {
    let (x, y) = at(s)
    beta_prefix(a, b, x, y) / (x * y)
  }
  let g = (u : Double) => {
    let s = @math.exp(u)
    let v = value(s)
    (@math.ln(v) - lt, s * density(s) / v)
  }
  let s0 = if upper { 1.0 - x0 } else { x0 }
  let s0 = if s0 > 0.0 && s0 < 1.0 { s0 } else { 0.5 }
  let mut lo = s0
  let mut hi = s0
  let mut rounds = 0
  while value(hi) < target && hi < 1.0 && rounds < 2100 {
    hi = if hi * 2.0 >= 1.0 { 1.0 } else { hi * 2.0 }
    rounds = rounds + 1
  }
  rounds = 0
  while value(lo) >= target && lo > 1.0e-300 && rounds < 2100 {
    lo = lo * 0.5
    rounds = rounds + 1
  }
  let r = @num.find_root_newton(
    g,
    @math.ln(lo),
    @math.ln(hi),
    @math.ln(s0),
    max_iterations=300,
  )
  let mut s = @math.exp(r.value)
  for _ in 0..<4 {
    let step = (value(s) - target) / density(s)
    if step.is_nan() || step.is_inf() {
      break
    }
    s = s - step
    if step.abs() <= 2.220446049250313e-16 * s.abs() {
      break
    }
  }
  if upper {
    1.0 - s
  } else {
    s
  }
}

///|
/// Owen's T(h, a) = (1/2π) ∫_0^a exp(-h²(1+t²)/2) / (1+t²) dt. Computed by
/// double-exponential quadrature for |a| <= 1, and through
/// T(h, a) + T(ah, 1/a) = (Φ(h)Φc(ah) + Φ(ah)Φc(h))/2 (h >= 0) beyond, the
/// right side written without cancellation.
pub fn owens_t(h : Double, a : Double) -> Double {
  if h.is_nan() || a.is_nan() {
    return 0.0 / 0.0
  }
  if a < 0.0 {
    return -owens_t(h, -a)
  }
  if a == 0.0 {
    return 0.0
  }
  let h = h.abs()
  if h == 0.0 {
    return @math.atan(a) / (2.0 * @math.PI)
  }
  if a <= 1.0 {
    let r = @num.integrate(
      t => {
        let u = 1.0 + t * t
        @math.exp(-0.5 * h * h * u) / u
      },
      0.0,
      a,
      rel_tol=1.0e-15,
    )
    r.value / (2.0 * @math.PI)
  } else {
    let ah = a * h
    let ph = ndtr(h)
    let qh = ndtr(-h)
    let pah = ndtr(ah)
    let qah = ndtr(-ah)
    0.5 * (ph * qah + pah * qh) - owens_t(ah, 1.0 / a)
  }
}

///|
/// The upper tail of Kolmogorov's limiting distribution,
/// P(K > x) = 2 Σ_{k>=1} (-1)^{k-1} exp(-2k²x²), with the theta-function
/// form of the lower tail for small x where that sum converges slowly.
pub fn kolmogorov_sf(x : Double) -> Double {
  if x.is_nan() {
    return x
  }
  if x <= 0.0 {
    return 1.0
  }
  if x < 1.0 {
    // P(K <= x) = √(2π)/x Σ exp(-(2k-1)²π²/(8x²))
    let c = @math.PI * @math.PI / (8.0 * x * x)
    let mut s = 0.0
    for k in 1..<=20 {
      let m = (2 * k - 1).to_double()
      let t = @math.exp(-m * m * c)
      s = s + t
      if t < 1.0e-18 * s {
        break
      }
    }
    1.0 - (2.0 * @math.PI).sqrt() / x * s
  } else {
    let mut s = 0.0
    let mut sign = 1.0
    for k in 1..<=100 {
      let kd = k.to_double()
      let t = @math.exp(-2.0 * kd * kd * x * x)
      s = s + sign * t
      sign = -sign
      if t < 1.0e-18 * s.abs() {
        break
      }
    }
    2.0 * s
  }
}

///|
/// The upper tail of the studentized range distribution (Tukey's HSD) for
/// k groups and df degrees of freedom: a double integral, not yet written.
pub fn studentized_range_sf(q : Double, k : Double, df : Double) -> Double {
  ignore((q, k, df))
  ...
}
