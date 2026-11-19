///|
/// The regularised incomplete gamma functions P(a, x), Q(a, x) and their
/// inverses.
///
/// The prefactor x^a e^{-x} / Γ(a) is the usual source of lost digits: for
/// a in the thousands its logarithm is a difference of numbers near
/// a ln a. For a >= 8 it is formed instead as
/// √(a/2π) exp(a·(ln(1+t) - t) - δ(a)) with t = (x - a)/a and δ Stirling's
/// remainder, where ln(1+t) - t has its own cancellation-free series.

///|
/// ln(1 + t) - t, accurate for small |t|.
fn log1pmx(t : Double) -> Double {
  if t.abs() < 0.5 {
    // -t²/2 + t³/3 - t⁴/4 + ...
    let mut sum = 0.0
    let mut tk = t
    let mut k = 2
    while k < 200 {
      tk = tk * -t
      let term = tk / k.to_double()
      sum = sum + term
      if term.abs() <= 1.0e-17 * sum.abs() {
        break
      }
      k = k + 1
    }
    sum
  } else {
    @math.ln_1p(t) - t
  }
}

///|
/// x^a e^{-x} / Γ(a).
fn gamma_prefix(a : Double, x : Double) -> Double {
  if x == 0.0 {
    return 0.0
  }
  if a >= 8.0 {
    let t = (x - a) / a
    // ln(1+t) - t; far from a, 1 + t = x/a and ln1p(t) would cancel
    let l = if t.abs() < 0.5 { log1pmx(t) } else { @math.ln(x / a) - t }
    (a / (2.0 * @math.PI)).sqrt() * @math.exp(a * l - stirling_delta(a))
  } else if x < 700.0 {
    @math.exp(-x) * @math.exp(a * @math.ln(x) - lgamma(a))
  } else {
    @math.exp(a * @math.ln(x) - x - lgamma(a))
  }
}

///|
/// P by its series: x^a e^{-x}/Γ(a+1) Σ x^n / ((a+1)...(a+n)).
fn igam_series(a : Double, x : Double) -> Double {
  let mut term = 1.0
  let mut sum = 1.0
  let mut n = 0
  while n < 100000 {
    n = n + 1
    term = term * x / (a + n.to_double())
    sum = sum + term
    if term <= 1.0e-17 * sum {
      break
    }
  }
  gamma_prefix(a, x) / a * sum
}

///|
/// Q by Legendre's continued fraction, for x >= a + 1.
fn igamc_cf(a : Double, x : Double) -> Double {
  let r = @num.continued_fraction(
    x + 1.0 - a,
    n => {
      let nd = n.to_double()
      (-nd * (nd - a), x + 2.0 * nd + 1.0 - a)
    },
    max_terms=100000,
  )
  gamma_prefix(a, x) / r.value
}

///|
/// Q for a < 1 and x < a + 1, where 1 - P would cancel:
/// Q = -expm1(a ln x - ln Γ(1+a)) - (x^a/Γ(a)) Σ_{n>=1} (-x)^n / (n! (a+n)).
fn igamc_small_a(a : Double, x : Double) -> Double {
  let lx = @math.ln(x)
  let head = -@math.expm1(a * lx - lgamma_1p(a))
  let xa_over_gamma = a * @math.exp(a * lx - lgamma_1p(a))
  let mut term = 1.0
  let mut sum = 0.0
  let mut n = 0
  while n < 1000 {
    n = n + 1
    term = term * -x / n.to_double()
    let t = term / (a + n.to_double())
    sum = sum + t
    if t.abs() <= 1.0e-17 * sum.abs() {
      break
    }
  }
  head - xa_over_gamma * sum
}

///|
/// (P(a, x), Q(a, x)), each computed directly where it is the small one.
pub fn gammainc_pair(a : Double, x : Double) -> (Double, Double) {
  if a.is_nan() || x.is_nan() || a <= 0.0 || x < 0.0 {
    return (0.0 / 0.0, 0.0 / 0.0)
  }
  if x == 0.0 {
    return (0.0, 1.0)
  }
  if x.is_inf() {
    return (1.0, 0.0)
  }
  if x < a + 1.0 {
    let p = igam_series(a, x)
    if a < 1.0 {
      (p, igamc_small_a(a, x))
    } else {
      (p, 1.0 - p)
    }
  } else {
    let q = igamc_cf(a, x)
    (1.0 - q, q)
  }
}

///|
/// The lower regularised incomplete gamma P(a, x), a > 0, x >= 0.
pub fn gammainc(a : Double, x : Double) -> Double {
  gammainc_pair(a, x).0
}

///|
/// The upper regularised incomplete gamma Q(a, x) = 1 - P(a, x).
pub fn gammaincc(a : Double, x : Double) -> Double {
  gammainc_pair(a, x).1
}

///|
/// Solve P(a, x) = target (lower) or Q(a, x) = target (upper) for x, with
/// Newton steps on the logarithm inside a bracket.
fn igam_invert(a : Double, target : Double, upper : Bool) -> Double {
  // the tail to solve in is the one given; its target is at most 1/2
  let p_lower = if upper { 1.0 - target } else { target }
  // start: Wilson–Hilferty for a >= 1, the small-x power law otherwise
  let z = ndtri(p_lower)
  let wh = a *
    @math.pow(1.0 - 1.0 / (9.0 * a) + z / (3.0 * a.sqrt()), 3.0)
  let power = @math.exp((@math.ln(p_lower) + lgamma(a + 1.0)) / a)
  let mut x0 = if a >= 1.0 && wh > 0.0 { wh } else { power }
  if !(x0 > 0.0) || x0.is_inf() {
    x0 = a
  }
  let lt = @math.ln(target)
  let value = (x : Double) => {
    let (p, q) = gammainc_pair(a, x)
    if upper {
      q
    } else {
      p
    }
  }
  // ln P (or ln Q) is nearly linear in u = ln x, so Newton steps in u
  // converge in a few iterations at any scale
  let g = (u : Double) => {
    let x = @math.exp(u)
    let v = value(x)
    let d = gamma_prefix(a, x) / v
    if upper {
      (@math.ln(v) - lt, -d)
    } else {
      (@math.ln(v) - lt, d)
    }
  }
  // bracket the root around the start, multiplicatively
  let mut lo = x0
  let mut hi = x0
  let below = (x : Double) => if upper { value(x) > target } else { value(x) < target }
  let mut rounds = 0
  while below(hi) && rounds < 2100 {
    hi = hi * 2.0
    rounds = rounds + 1
  }
  rounds = 0
  while !below(lo) && lo > 1.0e-300 && rounds < 2100 {
    lo = lo * 0.5
    rounds = rounds + 1
  }
  let r = @num.find_root_newton(g, @math.ln(lo), @math.ln(hi), @math.ln(x0), max_iterations=200)
  // exp(u) turns the rounding of u into a relative error |u|·eps in x;
  // finish with Newton steps on the value itself in x
  let mut x = @math.exp(r.value)
  for _ in 0..<4 {
    let v = value(x)
    let dv = gamma_prefix(a, x) / x
    let step = (v - target) / (if upper { -dv } else { dv })
    if step.is_nan() || step.is_inf() {
      break
    }
    x = x - step
    if step.abs() <= 2.220446049250313e-16 * x.abs() {
      break
    }
  }
  x
}

///|
/// x with P(a, x) = p.
pub fn gammaincinv(a : Double, p : Double) -> Double {
  if a.is_nan() || p.is_nan() || a <= 0.0 || p < 0.0 || p > 1.0 {
    return 0.0 / 0.0
  }
  if p == 0.0 {
    return 0.0
  }
  if p == 1.0 {
    return 1.0 / 0.0
  }
  if p <= 0.5 {
    igam_invert(a, p, false)
  } else {
    igam_invert(a, 1.0 - p, true)
  }
}

///|
/// x with Q(a, x) = q.
pub fn gammainccinv(a : Double, q : Double) -> Double {
  if a.is_nan() || q.is_nan() || a <= 0.0 || q < 0.0 || q > 1.0 {
    return 0.0 / 0.0
  }
  if q == 0.0 {
    return 1.0 / 0.0
  }
  if q == 1.0 {
    return 0.0
  }
  if q <= 0.5 {
    igam_invert(a, q, true)
  } else {
    igam_invert(a, 1.0 - q, false)
  }
}
