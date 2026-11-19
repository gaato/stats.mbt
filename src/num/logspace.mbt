///|
/// Arithmetic on log-probabilities. Tail probabilities below 1e-308 and
/// likelihoods of thousands of points only exist on the log scale, and the
/// naive formulas lose them: log(exp(a) + exp(b)) overflows, log(1 - exp(a))
/// cancels near a = 0. Each function here is the stable form.

///|
/// log(exp(a) + exp(b)).
pub fn log_add_exp(a : Double, b : Double) -> Double {
  if a == -1.0 / 0.0 {
    return b
  }
  if b == -1.0 / 0.0 {
    return a
  }
  let (hi, lo) = if a > b { (a, b) } else { (b, a) }
  hi + log1p(@math.exp(lo - hi))
}

///|
/// log(exp(a) - exp(b)) for a >= b; -inf when they are equal.
pub fn log_sub_exp(a : Double, b : Double) -> Double {
  if b == -1.0 / 0.0 {
    return a
  }
  a + log1m_exp(b - a)
}

///|
/// log(1 - exp(x)) for x <= 0, switching formula at -ln 2 (Mächler 2012):
/// log(-expm1(x)) near 0, log1p(-exp(x)) further out.
pub fn log1m_exp(x : Double) -> Double {
  if x > -0.6931471805599453 {
    @math.ln(-@math.expm1(x))
  } else {
    log1p(-@math.exp(x))
  }
}

///|
/// log(1 + exp(x)), the softplus, without overflow for large x.
pub fn log1p_exp(x : Double) -> Double {
  if x > 33.3 {
    x
  } else if x > -37.0 {
    log1p(@math.exp(x))
  } else {
    @math.exp(x)
  }
}

///|
/// log Σ exp(xᵢ), shifted by the maximum. -inf for an empty input.
pub fn log_sum_exp(xs : ArrayView[Double]) -> Double {
  let mut m = -1.0 / 0.0
  for x in xs {
    if x > m {
      m = x
    }
  }
  if m == -1.0 / 0.0 || m == 1.0 / 0.0 {
    return m
  }
  let mut s = 0.0
  for x in xs {
    s = s + @math.exp(x - m)
  }
  m + @math.ln(s)
}

///|
/// log(1 + x), accurate near 0 (core's `@math` has `expm1` but this
/// package keeps its own `log1p` so the whole log-space layer is in one
/// place; replace with the core function once it is confirmed accurate).
pub fn log1p(x : Double) -> Double {
  let u = 1.0 + x
  if u == 1.0 {
    x
  } else {
    @math.ln(u) * x / (u - 1.0)
  }
}
