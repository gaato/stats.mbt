///|
/// Closed forms for the catalogue families that need no special functions.
/// Every formula is written for its tails: upper tails come from their own
/// expression (never 1 - cdf), and expm1, ln_1p and atan(1/x) replace the
/// forms that cancel. Checked against SciPy by catalog_test.mbt.

///|
fn ln1p(x : Double) -> Double {
  @math.ln_1p(x)
}

///|
fn neg_inf() -> Double {
  -1.0 / 0.0
}

// ---------------------------------------------------------------------------
// Uniform, exponential, Laplace, logistic, Cauchy
// ---------------------------------------------------------------------------

///|
fn Uniform::std_pdf(_self : Uniform, x : Double) -> Double {
  if x < 0.0 || x > 1.0 {
    0.0
  } else {
    1.0
  }
}

///|
fn Uniform::std_cdf(_self : Uniform, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    x
  }
}

///|
fn Uniform::std_sf(_self : Uniform, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else {
    1.0 - x
  }
}

///|
fn Uniform::std_quantile(_self : Uniform, p : Double) -> Double {
  p
}

///|
fn Uniform::std_isf(_self : Uniform, q : Double) -> Double {
  1.0 - q
}

///|
fn Exponential::std_pdf(_self : Exponential, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(-x)
  }
}

///|
fn Exponential::std_logpdf(_self : Exponential, x : Double) -> Double {
  if x < 0.0 {
    neg_inf()
  } else {
    -x
  }
}

///|
fn Exponential::std_cdf(_self : Exponential, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-x)
  }
}

///|
fn Exponential::std_sf(_self : Exponential, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-x)
  }
}

///|
fn Exponential::std_quantile(_self : Exponential, p : Double) -> Double {
  -ln1p(-p)
}

///|
fn Exponential::std_isf(_self : Exponential, q : Double) -> Double {
  -@math.ln(q)
}

///|
fn Laplace::std_pdf(_self : Laplace, x : Double) -> Double {
  0.5 * @math.exp(-x.abs())
}

///|
fn Laplace::std_logpdf(_self : Laplace, x : Double) -> Double {
  -@math.ln(2.0) - x.abs()
}

///|
fn Laplace::std_cdf(_self : Laplace, x : Double) -> Double {
  if x < 0.0 {
    0.5 * @math.exp(x)
  } else {
    1.0 - 0.5 * @math.exp(-x)
  }
}

///|
fn Laplace::std_sf(_self : Laplace, x : Double) -> Double {
  if x > 0.0 {
    0.5 * @math.exp(-x)
  } else {
    1.0 - 0.5 * @math.exp(x)
  }
}

///|
fn Laplace::std_quantile(_self : Laplace, p : Double) -> Double {
  if p <= 0.5 {
    @math.ln(2.0 * p)
  } else {
    -@math.ln(2.0 * (1.0 - p))
  }
}

///|
fn Laplace::std_isf(self : Laplace, q : Double) -> Double {
  -self.std_quantile(q)
}

///|
fn Logistic::std_pdf(self : Logistic, x : Double) -> Double {
  @math.exp(self.std_logpdf(x))
}

///|
fn Logistic::std_logpdf(_self : Logistic, x : Double) -> Double {
  let a = x.abs()
  -a - 2.0 * ln1p(@math.exp(-a))
}

///|
fn Logistic::std_cdf(_self : Logistic, x : Double) -> Double {
  if x >= 0.0 {
    1.0 / (1.0 + @math.exp(-x))
  } else {
    let e = @math.exp(x)
    e / (1.0 + e)
  }
}

///|
fn Logistic::std_sf(self : Logistic, x : Double) -> Double {
  self.std_cdf(-x)
}

///|
fn Logistic::std_quantile(_self : Logistic, p : Double) -> Double {
  @math.ln(p) - ln1p(-p)
}

///|
fn Logistic::std_isf(self : Logistic, q : Double) -> Double {
  -self.std_quantile(q)
}

///|
fn Cauchy::std_pdf(_self : Cauchy, x : Double) -> Double {
  1.0 / (@math.PI * (1.0 + x * x))
}

///|
/// 1/2 + atan(x)/π, written as atan(-1/x)/π for x < 0 so the lower tail
/// keeps its digits.
fn Cauchy::std_cdf(_self : Cauchy, x : Double) -> Double {
  if x < 0.0 {
    @math.atan(-1.0 / x) / @math.PI
  } else if x == 0.0 {
    0.5
  } else {
    1.0 - @math.atan(1.0 / x) / @math.PI
  }
}

///|
fn Cauchy::std_sf(self : Cauchy, x : Double) -> Double {
  self.std_cdf(-x)
}

///|
fn Cauchy::std_quantile(_self : Cauchy, p : Double) -> Double {
  if p < 0.5 {
    -1.0 / @math.tan(@math.PI * p)
  } else if p == 0.5 {
    0.0
  } else {
    1.0 / @math.tan(@math.PI * (1.0 - p))
  }
}

///|
fn Cauchy::std_isf(self : Cauchy, q : Double) -> Double {
  -self.std_quantile(q)
}

// ---------------------------------------------------------------------------
// Gumbel, Weibull, Pareto, Lomax, Rayleigh
// ---------------------------------------------------------------------------

///|
fn GumbelRight::std_pdf(self : GumbelRight, x : Double) -> Double {
  @math.exp(self.std_logpdf(x))
}

///|
fn GumbelRight::std_logpdf(_self : GumbelRight, x : Double) -> Double {
  -x - @math.exp(-x)
}

///|
fn GumbelRight::std_cdf(_self : GumbelRight, x : Double) -> Double {
  @math.exp(-@math.exp(-x))
}

///|
fn GumbelRight::std_sf(_self : GumbelRight, x : Double) -> Double {
  -@math.expm1(-@math.exp(-x))
}

///|
fn GumbelRight::std_quantile(_self : GumbelRight, p : Double) -> Double {
  -@math.ln(-@math.ln(p))
}

///|
fn GumbelRight::std_isf(_self : GumbelRight, q : Double) -> Double {
  -@math.ln(-ln1p(-q))
}

///|
fn GumbelLeft::std_pdf(self : GumbelLeft, x : Double) -> Double {
  @math.exp(self.std_logpdf(x))
}

///|
fn GumbelLeft::std_logpdf(_self : GumbelLeft, x : Double) -> Double {
  x - @math.exp(x)
}

///|
fn GumbelLeft::std_cdf(_self : GumbelLeft, x : Double) -> Double {
  -@math.expm1(-@math.exp(x))
}

///|
fn GumbelLeft::std_sf(_self : GumbelLeft, x : Double) -> Double {
  @math.exp(-@math.exp(x))
}

///|
fn GumbelLeft::std_quantile(_self : GumbelLeft, p : Double) -> Double {
  @math.ln(-ln1p(-p))
}

///|
fn GumbelLeft::std_isf(_self : GumbelLeft, q : Double) -> Double {
  @math.ln(-@math.ln(q))
}

///|
fn Weibull::std_pdf(self : Weibull, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Weibull::std_logpdf(self : Weibull, x : Double) -> Double {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) + (self.c - 1.0) * @math.ln(x) - @math.pow(x, self.c)
  }
}

///|
fn Weibull::std_cdf(self : Weibull, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-@math.pow(x, self.c))
  }
}

///|
fn Weibull::std_sf(self : Weibull, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-@math.pow(x, self.c))
  }
}

///|
fn Weibull::std_quantile(self : Weibull, p : Double) -> Double {
  @math.pow(-ln1p(-p), 1.0 / self.c)
}

///|
fn Weibull::std_isf(self : Weibull, q : Double) -> Double {
  @math.pow(-@math.ln(q), 1.0 / self.c)
}

///|
fn WeibullMax::std_pdf(self : WeibullMax, x : Double) -> Double {
  if x > 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn WeibullMax::std_logpdf(self : WeibullMax, x : Double) -> Double {
  if x > 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) + (self.c - 1.0) * @math.ln(-x) - @math.pow(-x, self.c)
  }
}

///|
fn WeibullMax::std_cdf(self : WeibullMax, x : Double) -> Double {
  if x >= 0.0 {
    1.0
  } else {
    @math.exp(-@math.pow(-x, self.c))
  }
}

///|
fn WeibullMax::std_sf(self : WeibullMax, x : Double) -> Double {
  if x >= 0.0 {
    0.0
  } else {
    -@math.expm1(-@math.pow(-x, self.c))
  }
}

///|
fn WeibullMax::std_quantile(self : WeibullMax, p : Double) -> Double {
  -@math.pow(-@math.ln(p), 1.0 / self.c)
}

///|
fn WeibullMax::std_isf(self : WeibullMax, q : Double) -> Double {
  -@math.pow(-ln1p(-q), 1.0 / self.c)
}

///|
fn Pareto::std_pdf(self : Pareto, x : Double) -> Double {
  if x < 1.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Pareto::std_logpdf(self : Pareto, x : Double) -> Double {
  if x < 1.0 {
    neg_inf()
  } else {
    @math.ln(self.b) - (self.b + 1.0) * @math.ln(x)
  }
}

///|
fn Pareto::std_cdf(self : Pareto, x : Double) -> Double {
  if x <= 1.0 {
    0.0
  } else {
    -@math.expm1(-self.b * @math.ln(x))
  }
}

///|
fn Pareto::std_sf(self : Pareto, x : Double) -> Double {
  if x <= 1.0 {
    1.0
  } else {
    @math.exp(-self.b * @math.ln(x))
  }
}

///|
fn Pareto::std_quantile(self : Pareto, p : Double) -> Double {
  @math.exp(-ln1p(-p) / self.b)
}

///|
fn Pareto::std_isf(self : Pareto, q : Double) -> Double {
  @math.exp(-@math.ln(q) / self.b)
}

///|
fn Lomax::std_pdf(self : Lomax, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Lomax::std_logpdf(self : Lomax, x : Double) -> Double {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) - (self.c + 1.0) * ln1p(x)
  }
}

///|
fn Lomax::std_cdf(self : Lomax, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-self.c * ln1p(x))
  }
}

///|
fn Lomax::std_sf(self : Lomax, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-self.c * ln1p(x))
  }
}

///|
fn Lomax::std_quantile(self : Lomax, p : Double) -> Double {
  @math.expm1(-ln1p(-p) / self.c)
}

///|
fn Lomax::std_isf(self : Lomax, q : Double) -> Double {
  @math.expm1(-@math.ln(q) / self.c)
}

///|
fn Rayleigh::std_pdf(_self : Rayleigh, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    x * @math.exp(-0.5 * x * x)
  }
}

///|
fn Rayleigh::std_logpdf(_self : Rayleigh, x : Double) -> Double {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(x) - 0.5 * x * x
  }
}

///|
fn Rayleigh::std_cdf(_self : Rayleigh, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-0.5 * x * x)
  }
}

///|
fn Rayleigh::std_sf(_self : Rayleigh, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-0.5 * x * x)
  }
}

///|
fn Rayleigh::std_quantile(_self : Rayleigh, p : Double) -> Double {
  (-2.0 * ln1p(-p)).sqrt()
}

///|
fn Rayleigh::std_isf(_self : Rayleigh, q : Double) -> Double {
  (-2.0 * @math.ln(q)).sqrt()
}

// ---------------------------------------------------------------------------
// Bounded: triangular, arcsine, power law, Bradford, log-uniform, anglit,
// truncated exponential, truncated Pareto
// ---------------------------------------------------------------------------

///|
fn Triangular::std_pdf(self : Triangular, x : Double) -> Double {
  if x < 0.0 || x > 1.0 {
    0.0
  } else if x < self.c {
    2.0 * x / self.c
  } else {
    2.0 * (1.0 - x) / (1.0 - self.c)
  }
}

///|
fn Triangular::std_cdf(self : Triangular, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else if x < self.c {
    x * x / self.c
  } else {
    1.0 - (1.0 - x) * (1.0 - x) / (1.0 - self.c)
  }
}

///|
fn Triangular::std_sf(self : Triangular, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else if x < self.c {
    1.0 - x * x / self.c
  } else {
    (1.0 - x) * (1.0 - x) / (1.0 - self.c)
  }
}

///|
fn Triangular::std_quantile(self : Triangular, p : Double) -> Double {
  if p < self.c {
    (self.c * p).sqrt()
  } else {
    1.0 - ((1.0 - self.c) * (1.0 - p)).sqrt()
  }
}

///|
fn Triangular::std_isf(self : Triangular, q : Double) -> Double {
  if q > 1.0 - self.c {
    (self.c * (1.0 - q)).sqrt()
  } else {
    1.0 - ((1.0 - self.c) * q).sqrt()
  }
}

///|
fn Arcsine::std_pdf(_self : Arcsine, x : Double) -> Double {
  if x <= 0.0 || x >= 1.0 {
    0.0
  } else {
    1.0 / (@math.PI * (x * (1.0 - x)).sqrt())
  }
}

///|
fn Arcsine::std_cdf(_self : Arcsine, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else if x <= 0.5 {
    2.0 / @math.PI * @math.asin(x.sqrt())
  } else {
    2.0 / @math.PI * @math.acos((1.0 - x).sqrt())
  }
}

///|
/// (2/π) acos(√x), switching to asin(√(1 - x)) above 1/2 where acos loses
/// digits; each branch keeps its argument away from where it is
/// ill-conditioned.
fn Arcsine::std_sf(_self : Arcsine, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else if x < 0.5 {
    2.0 / @math.PI * @math.acos(x.sqrt())
  } else {
    2.0 / @math.PI * @math.asin((1.0 - x).sqrt())
  }
}

///|
fn Arcsine::std_quantile(_self : Arcsine, p : Double) -> Double {
  let s = @math.sin(0.5 * @math.PI * p)
  s * s
}

///|
fn Arcsine::std_isf(_self : Arcsine, q : Double) -> Double {
  let c = @math.cos(0.5 * @math.PI * q)
  c * c
}

///|
fn PowerLaw::std_pdf(self : PowerLaw, x : Double) -> Double {
  if x < 0.0 || x > 1.0 {
    0.0
  } else {
    self.a * @math.pow(x, self.a - 1.0)
  }
}

///|
fn PowerLaw::std_cdf(self : PowerLaw, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    @math.pow(x, self.a)
  }
}

///|
fn PowerLaw::std_sf(self : PowerLaw, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else {
    -@math.expm1(self.a * @math.ln(x))
  }
}

///|
fn PowerLaw::std_quantile(self : PowerLaw, p : Double) -> Double {
  @math.pow(p, 1.0 / self.a)
}

///|
fn PowerLaw::std_isf(self : PowerLaw, q : Double) -> Double {
  @math.exp(ln1p(-q) / self.a)
}

///|
fn Bradford::std_pdf(self : Bradford, x : Double) -> Double {
  if x < 0.0 || x > 1.0 {
    0.0
  } else {
    self.c / (ln1p(self.c) * (1.0 + self.c * x))
  }
}

///|
fn Bradford::std_cdf(self : Bradford, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    ln1p(self.c * x) / ln1p(self.c)
  }
}

///|
fn Bradford::std_sf(self : Bradford, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else {
    // ln((1 + c)/(1 + cx)) = ln1p(c(1 - x)/(1 + cx))
    ln1p(self.c * (1.0 - x) / (1.0 + self.c * x)) / ln1p(self.c)
  }
}

///|
fn Bradford::std_quantile(self : Bradford, p : Double) -> Double {
  @math.expm1(p * ln1p(self.c)) / self.c
}

///|
fn Bradford::std_isf(self : Bradford, q : Double) -> Double {
  @math.expm1((1.0 - q) * ln1p(self.c)) / self.c
}

///|
fn LogUniform::std_pdf(self : LogUniform, x : Double) -> Double {
  if x < self.a || x > self.b {
    0.0
  } else {
    1.0 / (x * (@math.ln(self.b) - @math.ln(self.a)))
  }
}

///|
fn LogUniform::std_cdf(self : LogUniform, x : Double) -> Double {
  if x <= self.a {
    0.0
  } else if x >= self.b {
    1.0
  } else {
    // ln(x/a) = ln1p((x - a)/a); x - a is exact near a
    ln1p((x - self.a) / self.a) / (@math.ln(self.b) - @math.ln(self.a))
  }
}

///|
fn LogUniform::std_sf(self : LogUniform, x : Double) -> Double {
  if x <= self.a {
    1.0
  } else if x >= self.b {
    0.0
  } else {
    ln1p((self.b - x) / x) / (@math.ln(self.b) - @math.ln(self.a))
  }
}

///|
fn LogUniform::std_quantile(self : LogUniform, p : Double) -> Double {
  @math.exp(@math.ln(self.a) + p * (@math.ln(self.b) - @math.ln(self.a)))
}

///|
fn LogUniform::std_isf(self : LogUniform, q : Double) -> Double {
  @math.exp(@math.ln(self.b) - q * (@math.ln(self.b) - @math.ln(self.a)))
}

///|
fn Anglit::std_pdf(_self : Anglit, x : Double) -> Double {
  if x < -0.25 * @math.PI || x > 0.25 * @math.PI {
    0.0
  } else {
    @math.cos(2.0 * x)
  }
}

///|
fn Anglit::std_cdf(_self : Anglit, x : Double) -> Double {
  if x <= -0.25 * @math.PI {
    0.0
  } else if x >= 0.25 * @math.PI {
    1.0
  } else {
    let s = @math.sin(x + 0.25 * @math.PI)
    s * s
  }
}

///|
fn Anglit::std_sf(_self : Anglit, x : Double) -> Double {
  if x <= -0.25 * @math.PI {
    1.0
  } else if x >= 0.25 * @math.PI {
    0.0
  } else {
    let c = @math.cos(x + 0.25 * @math.PI)
    c * c
  }
}

///|
fn Anglit::std_quantile(_self : Anglit, p : Double) -> Double {
  @math.asin(p.sqrt()) - 0.25 * @math.PI
}

///|
fn Anglit::std_isf(_self : Anglit, q : Double) -> Double {
  @math.acos(q.sqrt()) - 0.25 * @math.PI
}

///|
fn TruncatedExponential::std_pdf(
  self : TruncatedExponential,
  x : Double,
) -> Double {
  if x < 0.0 || x > self.b {
    0.0
  } else {
    @math.exp(-x) / -@math.expm1(-self.b)
  }
}

///|
fn TruncatedExponential::std_cdf(
  self : TruncatedExponential,
  x : Double,
) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= self.b {
    1.0
  } else {
    @math.expm1(-x) / @math.expm1(-self.b)
  }
}

///|
/// (e^-x - e^-b)/(1 - e^-b) = e^-x (1 - e^-(b-x)) / (1 - e^-b).
fn TruncatedExponential::std_sf(
  self : TruncatedExponential,
  x : Double,
) -> Double {
  if x <= 0.0 {
    1.0
  } else if x >= self.b {
    0.0
  } else {
    @math.exp(-x) * @math.expm1(x - self.b) / @math.expm1(-self.b)
  }
}

///|
fn TruncatedExponential::std_quantile(
  self : TruncatedExponential,
  p : Double,
) -> Double {
  -ln1p(p * @math.expm1(-self.b))
}

///|
fn TruncatedExponential::std_isf(
  self : TruncatedExponential,
  q : Double,
) -> Double {
  -@math.ln(@math.exp(-self.b) - q * @math.expm1(-self.b))
}

///|
fn TruncatedPareto::mass(self : TruncatedPareto) -> Double {
  -@math.expm1(-self.b * @math.ln(self.c))
}

///|
fn TruncatedPareto::std_pdf(self : TruncatedPareto, x : Double) -> Double {
  if x < 1.0 || x > self.c {
    0.0
  } else {
    self.b * @math.exp(-(self.b + 1.0) * @math.ln(x)) / self.mass()
  }
}

///|
fn TruncatedPareto::std_cdf(self : TruncatedPareto, x : Double) -> Double {
  if x <= 1.0 {
    0.0
  } else if x >= self.c {
    1.0
  } else {
    -@math.expm1(-self.b * @math.ln(x)) / self.mass()
  }
}

///|
fn TruncatedPareto::std_sf(self : TruncatedPareto, x : Double) -> Double {
  if x <= 1.0 {
    1.0
  } else if x >= self.c {
    0.0
  } else {
    @math.exp(-self.b * @math.ln(x)) *
    -@math.expm1(-self.b * (@math.ln(self.c) - @math.ln(x))) /
    self.mass()
  }
}

///|
fn TruncatedPareto::std_quantile(self : TruncatedPareto, p : Double) -> Double {
  @math.exp(-ln1p(-p * self.mass()) / self.b)
}

///|
fn TruncatedPareto::std_isf(self : TruncatedPareto, q : Double) -> Double {
  @math.exp(
    -@math.ln(@math.exp(-self.b * @math.ln(self.c)) + q * self.mass()) / self.b,
  )
}

// ---------------------------------------------------------------------------
// Half and folded forms, hyperbolic secant
// ---------------------------------------------------------------------------

///|
fn HalfCauchy::std_pdf(_self : HalfCauchy, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    2.0 / (@math.PI * (1.0 + x * x))
  }
}

///|
fn HalfCauchy::std_cdf(_self : HalfCauchy, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    2.0 / @math.PI * @math.atan(x)
  }
}

///|
fn HalfCauchy::std_sf(_self : HalfCauchy, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    2.0 / @math.PI * @math.atan(1.0 / x)
  }
}

///|
/// tan(πp/2) is ill-conditioned near p = 1; use 1/tan(π(1 - p)/2) there.
fn HalfCauchy::std_quantile(_self : HalfCauchy, p : Double) -> Double {
  if p <= 0.5 {
    @math.tan(0.5 * @math.PI * p)
  } else {
    1.0 / @math.tan(0.5 * @math.PI * (1.0 - p))
  }
}

///|
fn HalfCauchy::std_isf(_self : HalfCauchy, q : Double) -> Double {
  1.0 / @math.tan(0.5 * @math.PI * q)
}

///|
fn HalfLogistic::std_pdf(self : HalfLogistic, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn HalfLogistic::std_logpdf(_self : HalfLogistic, x : Double) -> Double {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(2.0) - x - 2.0 * ln1p(@math.exp(-x))
  }
}

///|
fn HalfLogistic::std_cdf(_self : HalfLogistic, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.tanh(0.5 * x)
  }
}

///|
fn HalfLogistic::std_sf(_self : HalfLogistic, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    let e = @math.exp(-x)
    2.0 * e / (1.0 + e)
  }
}

///|
fn HalfLogistic::std_quantile(_self : HalfLogistic, p : Double) -> Double {
  ln1p(p) - ln1p(-p)
}

///|
/// ln((2 - q)/q) = ln1p(2(1 - q)/q), exact when q is near 1.
fn HalfLogistic::std_isf(_self : HalfLogistic, q : Double) -> Double {
  ln1p(2.0 * (1.0 - q) / q)
}

///|
fn HyperbolicSecant::std_pdf(_self : HyperbolicSecant, x : Double) -> Double {
  1.0 / (@math.PI * @math.cosh(x))
}

///|
fn HyperbolicSecant::std_cdf(_self : HyperbolicSecant, x : Double) -> Double {
  2.0 / @math.PI * @math.atan(@math.exp(x))
}

///|
fn HyperbolicSecant::std_sf(_self : HyperbolicSecant, x : Double) -> Double {
  2.0 / @math.PI * @math.atan(@math.exp(-x))
}

///|
fn HyperbolicSecant::std_quantile(
  _self : HyperbolicSecant,
  p : Double,
) -> Double {
  if p <= 0.5 {
    @math.ln(@math.tan(0.5 * @math.PI * p))
  } else {
    -@math.ln(@math.tan(0.5 * @math.PI * (1.0 - p)))
  }
}

///|
fn HyperbolicSecant::std_isf(self : HyperbolicSecant, q : Double) -> Double {
  -self.std_quantile(q)
}

// ---------------------------------------------------------------------------
// Generalised logistic, Pareto, extreme value
// ---------------------------------------------------------------------------

///|
fn GeneralizedLogistic::std_pdf(
  self : GeneralizedLogistic,
  x : Double,
) -> Double {
  @math.exp(self.std_logpdf(x))
}

///|
fn GeneralizedLogistic::std_logpdf(
  self : GeneralizedLogistic,
  x : Double,
) -> Double {
  @math.ln(self.c) - x - (self.c + 1.0) * @num.log1p_exp(-x)
}

///|
fn GeneralizedLogistic::std_cdf(
  self : GeneralizedLogistic,
  x : Double,
) -> Double {
  @math.exp(-self.c * @num.log1p_exp(-x))
}

///|
fn GeneralizedLogistic::std_sf(
  self : GeneralizedLogistic,
  x : Double,
) -> Double {
  -@math.expm1(-self.c * @num.log1p_exp(-x))
}

///|
fn GeneralizedLogistic::std_quantile(
  self : GeneralizedLogistic,
  p : Double,
) -> Double {
  -@math.ln(@math.expm1(-@math.ln(p) / self.c))
}

///|
fn GeneralizedLogistic::std_isf(
  self : GeneralizedLogistic,
  q : Double,
) -> Double {
  -@math.ln(@math.expm1(-ln1p(-q) / self.c))
}

///|
/// log S(x) = -ln(1 + cx)/c, or -x for c = 0; NaN-free inside the support.
fn GeneralizedPareto::log_sf(self : GeneralizedPareto, x : Double) -> Double {
  if self.c == 0.0 {
    -x
  } else {
    -ln1p(self.c * x) / self.c
  }
}

///|
fn GeneralizedPareto::upper(self : GeneralizedPareto) -> Double {
  if self.c < 0.0 {
    -1.0 / self.c
  } else {
    1.0 / 0.0
  }
}

///|
fn GeneralizedPareto::std_pdf(self : GeneralizedPareto, x : Double) -> Double {
  if x < 0.0 || x > self.upper() {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn GeneralizedPareto::std_logpdf(
  self : GeneralizedPareto,
  x : Double,
) -> Double {
  if x < 0.0 || x > self.upper() {
    neg_inf()
  } else if self.c == 0.0 {
    -x
  } else {
    -(1.0 + 1.0 / self.c) * ln1p(self.c * x)
  }
}

///|
fn GeneralizedPareto::std_cdf(self : GeneralizedPareto, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= self.upper() {
    1.0
  } else {
    -@math.expm1(self.log_sf(x))
  }
}

///|
fn GeneralizedPareto::std_sf(self : GeneralizedPareto, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else if x >= self.upper() {
    0.0
  } else {
    @math.exp(self.log_sf(x))
  }
}

///|
fn GeneralizedPareto::std_quantile(
  self : GeneralizedPareto,
  p : Double,
) -> Double {
  if self.c == 0.0 {
    -ln1p(-p)
  } else {
    @math.expm1(-self.c * ln1p(-p)) / self.c
  }
}

///|
fn GeneralizedPareto::std_isf(self : GeneralizedPareto, q : Double) -> Double {
  if self.c == 0.0 {
    -@math.ln(q)
  } else {
    @math.expm1(-self.c * @math.ln(q)) / self.c
  }
}

///|
/// t(x) = (1 - cx)^(1/c), or e^-x for c = 0; F = exp(-t).
fn GeneralizedExtremeValue::t(
  self : GeneralizedExtremeValue,
  x : Double,
) -> Double {
  if self.c == 0.0 {
    @math.exp(-x)
  } else {
    @math.exp(ln1p(-self.c * x) / self.c)
  }
}

///|
fn GeneralizedExtremeValue::inside(
  self : GeneralizedExtremeValue,
  x : Double,
) -> Bool {
  self.c == 0.0 || self.c * x < 1.0
}

///|
fn GeneralizedExtremeValue::std_pdf(
  self : GeneralizedExtremeValue,
  x : Double,
) -> Double {
  if !self.inside(x) {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn GeneralizedExtremeValue::std_logpdf(
  self : GeneralizedExtremeValue,
  x : Double,
) -> Double {
  if !self.inside(x) {
    neg_inf()
  } else if self.c == 0.0 {
    -x - @math.exp(-x)
  } else {
    -self.t(x) + (1.0 / self.c - 1.0) * ln1p(-self.c * x)
  }
}

///|
fn GeneralizedExtremeValue::std_cdf(
  self : GeneralizedExtremeValue,
  x : Double,
) -> Double {
  if !self.inside(x) {
    if self.c > 0.0 {
      1.0
    } else {
      0.0
    }
  } else {
    @math.exp(-self.t(x))
  }
}

///|
fn GeneralizedExtremeValue::std_sf(
  self : GeneralizedExtremeValue,
  x : Double,
) -> Double {
  if !self.inside(x) {
    if self.c > 0.0 {
      0.0
    } else {
      1.0
    }
  } else {
    -@math.expm1(-self.t(x))
  }
}

///|
fn GeneralizedExtremeValue::std_quantile(
  self : GeneralizedExtremeValue,
  p : Double,
) -> Double {
  let l = -@math.ln(p)
  if self.c == 0.0 {
    -@math.ln(l)
  } else {
    -@math.expm1(self.c * @math.ln(l)) / self.c
  }
}

///|
fn GeneralizedExtremeValue::std_isf(
  self : GeneralizedExtremeValue,
  q : Double,
) -> Double {
  let l = -ln1p(-q)
  if self.c == 0.0 {
    -@math.ln(l)
  } else {
    -@math.expm1(self.c * @math.ln(l)) / self.c
  }
}

// ---------------------------------------------------------------------------
// Burr family
// ---------------------------------------------------------------------------

///|
fn burr3_log_cdf(c : Double, d : Double, x : Double) -> Double {
  -d * ln1p(@math.pow(x, -c))
}

///|
fn Burr::std_pdf(self : Burr, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Burr::std_logpdf(self : Burr, x : Double) -> Double {
  if x <= 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) +
    @math.ln(self.d) -
    (self.c + 1.0) * @math.ln(x) -
    (self.d + 1.0) * ln1p(@math.pow(x, -self.c))
  }
}

///|
fn Burr::std_cdf(self : Burr, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(burr3_log_cdf(self.c, self.d, x))
  }
}

///|
fn Burr::std_sf(self : Burr, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    -@math.expm1(burr3_log_cdf(self.c, self.d, x))
  }
}

///|
fn Burr::std_quantile(self : Burr, p : Double) -> Double {
  @math.pow(@math.expm1(-@math.ln(p) / self.d), -1.0 / self.c)
}

///|
fn Burr::std_isf(self : Burr, q : Double) -> Double {
  @math.pow(@math.expm1(-ln1p(-q) / self.d), -1.0 / self.c)
}

///|
fn Fisk::std_pdf(self : Fisk, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Fisk::std_logpdf(self : Fisk, x : Double) -> Double {
  if x <= 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) -
    (self.c + 1.0) * @math.ln(x) -
    2.0 * ln1p(@math.pow(x, -self.c))
  }
}

///|
fn Fisk::std_cdf(self : Fisk, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(burr3_log_cdf(self.c, 1.0, x))
  }
}

///|
fn Fisk::std_sf(self : Fisk, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    -@math.expm1(burr3_log_cdf(self.c, 1.0, x))
  }
}

///|
fn Fisk::std_quantile(self : Fisk, p : Double) -> Double {
  @math.pow(@math.expm1(-@math.ln(p)), -1.0 / self.c)
}

///|
fn Fisk::std_isf(self : Fisk, q : Double) -> Double {
  @math.pow(@math.expm1(-ln1p(-q)), -1.0 / self.c)
}

///|
fn Burr12::std_pdf(self : Burr12, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Burr12::std_logpdf(self : Burr12, x : Double) -> Double {
  if x <= 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) +
    @math.ln(self.d) +
    (self.c - 1.0) * @math.ln(x) -
    (self.d + 1.0) * ln1p(@math.pow(x, self.c))
  }
}

///|
fn Burr12::std_cdf(self : Burr12, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-self.d * ln1p(@math.pow(x, self.c)))
  }
}

///|
fn Burr12::std_sf(self : Burr12, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-self.d * ln1p(@math.pow(x, self.c)))
  }
}

///|
fn Burr12::std_quantile(self : Burr12, p : Double) -> Double {
  @math.pow(@math.expm1(-ln1p(-p) / self.d), 1.0 / self.c)
}

///|
fn Burr12::std_isf(self : Burr12, q : Double) -> Double {
  @math.pow(@math.expm1(-@math.ln(q) / self.d), 1.0 / self.c)
}

///|
/// F = x (a + x^a)^(-1/a) = exp(-ln(1 + a x^-a)/a).
fn Kappa3::std_pdf(self : Kappa3, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    self.a * @math.pow(self.a + @math.pow(x, self.a), -1.0 / self.a - 1.0)
  }
}

///|
fn Kappa3::std_cdf(self : Kappa3, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(-ln1p(self.a * @math.pow(x, -self.a)) / self.a)
  }
}

///|
fn Kappa3::std_sf(self : Kappa3, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    -@math.expm1(-ln1p(self.a * @math.pow(x, -self.a)) / self.a)
  }
}

///|
fn Kappa3::std_quantile(self : Kappa3, p : Double) -> Double {
  @math.pow(self.a / @math.expm1(-self.a * @math.ln(p)), 1.0 / self.a)
}

///|
fn Kappa3::std_isf(self : Kappa3, q : Double) -> Double {
  @math.pow(self.a / @math.expm1(-self.a * ln1p(-q)), 1.0 / self.a)
}

// ---------------------------------------------------------------------------
// Exponential variants, Gompertz, double Weibull, log-Laplace,
// asymmetric Laplace, skew Cauchy
// ---------------------------------------------------------------------------

///|
fn ExponentialPower::std_pdf(self : ExponentialPower, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn ExponentialPower::std_logpdf(self : ExponentialPower, x : Double) -> Double {
  if x <= 0.0 {
    neg_inf()
  } else {
    let xb = @math.pow(x, self.b)
    @math.ln(self.b) + (self.b - 1.0) * @math.ln(x) + xb - @math.expm1(xb)
  }
}

///|
fn ExponentialPower::std_cdf(self : ExponentialPower, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-@math.expm1(@math.pow(x, self.b)))
  }
}

///|
fn ExponentialPower::std_sf(self : ExponentialPower, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-@math.expm1(@math.pow(x, self.b)))
  }
}

///|
fn ExponentialPower::std_quantile(
  self : ExponentialPower,
  p : Double,
) -> Double {
  @math.pow(ln1p(-ln1p(-p)), 1.0 / self.b)
}

///|
fn ExponentialPower::std_isf(self : ExponentialPower, q : Double) -> Double {
  @math.pow(ln1p(-@math.ln(q)), 1.0 / self.b)
}

///|
fn ExponentiatedWeibull::std_pdf(
  self : ExponentiatedWeibull,
  x : Double,
) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn ExponentiatedWeibull::std_logpdf(
  self : ExponentiatedWeibull,
  x : Double,
) -> Double {
  if x <= 0.0 {
    neg_inf()
  } else {
    let xc = @math.pow(x, self.c)
    @math.ln(self.a) +
    @math.ln(self.c) +
    (self.a - 1.0) * @math.ln(-@math.expm1(-xc)) -
    xc +
    (self.c - 1.0) * @math.ln(x)
  }
}

///|
fn ExponentiatedWeibull::std_cdf(
  self : ExponentiatedWeibull,
  x : Double,
) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(self.a * @math.ln(-@math.expm1(-@math.pow(x, self.c))))
  }
}

///|
fn ExponentiatedWeibull::std_sf(
  self : ExponentiatedWeibull,
  x : Double,
) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    -@math.expm1(self.a * @math.ln(-@math.expm1(-@math.pow(x, self.c))))
  }
}

///|
/// x^c = -ln(1 - p^(1/a)), with 1 - p^(1/a) = -expm1(ln(p)/a) so that
/// p^(1/a) near 1 does not cancel.
fn ExponentiatedWeibull::std_quantile(
  self : ExponentiatedWeibull,
  p : Double,
) -> Double {
  @math.pow(-@math.ln(-@math.expm1(@math.ln(p) / self.a)), 1.0 / self.c)
}

///|
fn ExponentiatedWeibull::std_isf(
  self : ExponentiatedWeibull,
  q : Double,
) -> Double {
  @math.pow(-@math.ln(-@math.expm1(ln1p(-q) / self.a)), 1.0 / self.c)
}

///|
fn Gompertz::std_pdf(self : Gompertz, x : Double) -> Double {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(self.std_logpdf(x))
  }
}

///|
fn Gompertz::std_logpdf(self : Gompertz, x : Double) -> Double {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) + x - self.c * @math.expm1(x)
  }
}

///|
fn Gompertz::std_cdf(self : Gompertz, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-self.c * @math.expm1(x))
  }
}

///|
fn Gompertz::std_sf(self : Gompertz, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-self.c * @math.expm1(x))
  }
}

///|
fn Gompertz::std_quantile(self : Gompertz, p : Double) -> Double {
  ln1p(-ln1p(-p) / self.c)
}

///|
fn Gompertz::std_isf(self : Gompertz, q : Double) -> Double {
  ln1p(-@math.ln(q) / self.c)
}

///|
fn DoubleWeibull::std_pdf(self : DoubleWeibull, x : Double) -> Double {
  let a = x.abs()
  0.5 * self.c * @math.pow(a, self.c - 1.0) * @math.exp(-@math.pow(a, self.c))
}

///|
fn DoubleWeibull::std_cdf(self : DoubleWeibull, x : Double) -> Double {
  if x <= 0.0 {
    0.5 * @math.exp(-@math.pow(-x, self.c))
  } else {
    1.0 - 0.5 * @math.exp(-@math.pow(x, self.c))
  }
}

///|
fn DoubleWeibull::std_sf(self : DoubleWeibull, x : Double) -> Double {
  self.std_cdf(-x)
}

///|
fn DoubleWeibull::std_quantile(self : DoubleWeibull, p : Double) -> Double {
  if p <= 0.5 {
    -@math.pow(-@math.ln(2.0 * p), 1.0 / self.c)
  } else {
    @math.pow(-@math.ln(2.0 * (1.0 - p)), 1.0 / self.c)
  }
}

///|
fn DoubleWeibull::std_isf(self : DoubleWeibull, q : Double) -> Double {
  -self.std_quantile(q)
}

///|
fn LogLaplace::std_pdf(self : LogLaplace, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x < 1.0 {
    0.5 * self.c * @math.pow(x, self.c - 1.0)
  } else {
    0.5 * self.c * @math.pow(x, -self.c - 1.0)
  }
}

///|
fn LogLaplace::std_cdf(self : LogLaplace, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x < 1.0 {
    0.5 * @math.pow(x, self.c)
  } else {
    1.0 - 0.5 * @math.pow(x, -self.c)
  }
}

///|
fn LogLaplace::std_sf(self : LogLaplace, x : Double) -> Double {
  if x <= 0.0 {
    1.0
  } else if x < 1.0 {
    1.0 - 0.5 * @math.pow(x, self.c)
  } else {
    0.5 * @math.pow(x, -self.c)
  }
}

///|
fn LogLaplace::std_quantile(self : LogLaplace, p : Double) -> Double {
  if p < 0.5 {
    @math.pow(2.0 * p, 1.0 / self.c)
  } else {
    @math.pow(2.0 * (1.0 - p), -1.0 / self.c)
  }
}

///|
fn LogLaplace::std_isf(self : LogLaplace, q : Double) -> Double {
  if q > 0.5 {
    @math.pow(2.0 * (1.0 - q), 1.0 / self.c)
  } else {
    @math.pow(2.0 * q, -1.0 / self.c)
  }
}

///|
fn AsymmetricLaplace::std_pdf(self : AsymmetricLaplace, x : Double) -> Double {
  let k = self.kappa
  let e = if x >= 0.0 { @math.exp(-x * k) } else { @math.exp(x / k) }
  e / (k + 1.0 / k)
}

///|
fn AsymmetricLaplace::std_cdf(self : AsymmetricLaplace, x : Double) -> Double {
  let k2 = self.kappa * self.kappa
  if x >= 0.0 {
    1.0 - @math.exp(-x * self.kappa) / (1.0 + k2)
  } else {
    k2 * @math.exp(x / self.kappa) / (1.0 + k2)
  }
}

///|
fn AsymmetricLaplace::std_sf(self : AsymmetricLaplace, x : Double) -> Double {
  let k2 = self.kappa * self.kappa
  if x >= 0.0 {
    @math.exp(-x * self.kappa) / (1.0 + k2)
  } else {
    1.0 - k2 * @math.exp(x / self.kappa) / (1.0 + k2)
  }
}

///|
fn AsymmetricLaplace::std_quantile(
  self : AsymmetricLaplace,
  p : Double,
) -> Double {
  let k = self.kappa
  let k2 = k * k
  if p <= k2 / (1.0 + k2) {
    k * @math.ln(p * (1.0 + k2) / k2)
  } else {
    -@math.ln((1.0 - p) * (1.0 + k2)) / k
  }
}

///|
fn AsymmetricLaplace::std_isf(self : AsymmetricLaplace, q : Double) -> Double {
  let k = self.kappa
  let k2 = k * k
  if q <= 1.0 / (1.0 + k2) {
    -@math.ln(q * (1.0 + k2)) / k
  } else {
    k * @math.ln((1.0 - q) * (1.0 + k2) / k2)
  }
}

///|
fn SkewCauchy::std_pdf(self : SkewCauchy, x : Double) -> Double {
  let s = if x < 0.0 { 1.0 - self.a } else { 1.0 + self.a }
  1.0 / (@math.PI * (x * x / (s * s) + 1.0))
}

///|
fn SkewCauchy::std_cdf(self : SkewCauchy, x : Double) -> Double {
  let lo = 1.0 - self.a
  let hi = 1.0 + self.a
  if x < 0.0 {
    lo / @math.PI * @math.atan(-lo / x)
  } else if x == 0.0 {
    0.5 * lo
  } else {
    1.0 - hi / @math.PI * @math.atan(hi / x)
  }
}

///|
fn SkewCauchy::std_sf(self : SkewCauchy, x : Double) -> Double {
  let lo = 1.0 - self.a
  let hi = 1.0 + self.a
  if x > 0.0 {
    hi / @math.PI * @math.atan(hi / x)
  } else if x == 0.0 {
    0.5 * hi
  } else {
    1.0 - lo / @math.PI * @math.atan(-lo / x)
  }
}

///|
fn SkewCauchy::std_quantile(self : SkewCauchy, p : Double) -> Double {
  let lo = 1.0 - self.a
  let hi = 1.0 + self.a
  if p < 0.5 * lo {
    -lo / @math.tan(@math.PI * p / lo)
  } else if p == 0.5 * lo {
    0.0
  } else {
    hi / @math.tan(@math.PI * (1.0 - p) / hi)
  }
}

///|
fn SkewCauchy::std_isf(self : SkewCauchy, q : Double) -> Double {
  let lo = 1.0 - self.a
  let hi = 1.0 + self.a
  if q < 0.5 * hi {
    hi / @math.tan(@math.PI * q / hi)
  } else if q == 0.5 * hi {
    0.0
  } else {
    -lo / @math.tan(@math.PI * (1.0 - q) / lo)
  }
}

// ---------------------------------------------------------------------------
// Log tails in closed form. Where a tail is itself an exponential its log
// is written down directly, so it stays finite after the tail underflows
// (an exponential upper tail at x = 900 is e^-900). The other tail keeps
// the trait's default, the log of the tail.
// ---------------------------------------------------------------------------

///|
fn Exponential::std_logsf(_self : Exponential, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -x
  }
}

///|
fn Laplace::std_logcdf(_self : Laplace, x : Double) -> Double {
  if x < 0.0 {
    x - @math.ln(2.0)
  } else {
    ln1p(-0.5 * @math.exp(-x))
  }
}

///|
fn Laplace::std_logsf(self : Laplace, x : Double) -> Double {
  self.std_logcdf(-x)
}

///|
/// ln(1/(1 + e^-x)).
fn Logistic::std_logcdf(_self : Logistic, x : Double) -> Double {
  -@num.log1p_exp(-x)
}

///|
fn Logistic::std_logsf(_self : Logistic, x : Double) -> Double {
  -@num.log1p_exp(x)
}

///|
fn GumbelRight::std_logcdf(_self : GumbelRight, x : Double) -> Double {
  -@math.exp(-x)
}

///|
fn GumbelLeft::std_logsf(_self : GumbelLeft, x : Double) -> Double {
  -@math.exp(x)
}

///|
fn Weibull::std_logsf(self : Weibull, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -@math.pow(x, self.c)
  }
}

///|
fn WeibullMax::std_logcdf(self : WeibullMax, x : Double) -> Double {
  if x >= 0.0 {
    0.0
  } else {
    -@math.pow(-x, self.c)
  }
}

///|
fn Pareto::std_logsf(self : Pareto, x : Double) -> Double {
  if x <= 1.0 {
    0.0
  } else {
    -self.b * @math.ln(x)
  }
}

///|
fn Lomax::std_logsf(self : Lomax, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -self.c * ln1p(x)
  }
}

///|
fn Rayleigh::std_logsf(_self : Rayleigh, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    -0.5 * x * x
  }
}

///|
fn PowerLaw::std_logcdf(self : PowerLaw, x : Double) -> Double {
  if x <= 0.0 {
    neg_inf()
  } else if x >= 1.0 {
    0.0
  } else {
    self.a * @math.ln(x)
  }
}

///|
/// ln(2 e^-x / (1 + e^-x)) = ln 2 - ln(1 + e^x).
fn HalfLogistic::std_logsf(_self : HalfLogistic, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    @math.ln(2.0) - @num.log1p_exp(x)
  }
}

///|
fn GeneralizedLogistic::std_logcdf(
  self : GeneralizedLogistic,
  x : Double,
) -> Double {
  -self.c * @num.log1p_exp(-x)
}

///|
fn GeneralizedPareto::std_logsf(self : GeneralizedPareto, x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x >= self.upper() {
    neg_inf()
  } else {
    self.log_sf(x)
  }
}

// ---------------------------------------------------------------------------
// Supports that depend on the shape parameters. A support that is the same
// for every member of a family is generated from SciPy's.
// ---------------------------------------------------------------------------

///|
fn GeneralizedPareto::std_support(self : GeneralizedPareto) -> (Double, Double) {
  (0.0, self.upper())
}

///|
/// Bounded above at 1/c for c > 0, below at 1/c for c < 0, the whole line
/// for c = 0.
fn GeneralizedExtremeValue::std_support(
  self : GeneralizedExtremeValue,
) -> (Double, Double) {
  if self.c > 0.0 {
    (neg_inf(), 1.0 / self.c)
  } else if self.c < 0.0 {
    (1.0 / self.c, 1.0 / 0.0)
  } else {
    (neg_inf(), 1.0 / 0.0)
  }
}

///|
fn LogUniform::std_support(self : LogUniform) -> (Double, Double) {
  (self.a, self.b)
}

///|
fn TruncatedExponential::std_support(
  self : TruncatedExponential,
) -> (Double, Double) {
  (0.0, self.b)
}

///|
fn TruncatedPareto::std_support(self : TruncatedPareto) -> (Double, Double) {
  (1.0, self.c)
}
