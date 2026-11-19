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
pub impl Density for Uniform with fn pdf(_self, x) {
  if x < 0.0 || x > 1.0 {
    0.0
  } else {
    1.0
  }
}

///|
pub impl Cdf for Uniform with fn cdf(_self, x) {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    x
  }
}

///|
pub impl Cdf for Uniform with fn sf(_self, x) {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else {
    1.0 - x
  }
}

///|
pub impl Quantile for Uniform with fn quantile(_self, p) {
  p
}

///|
pub impl Quantile for Uniform with fn isf(_self, q) {
  1.0 - q
}

///|
pub impl Density for Exponential with fn pdf(_self, x) {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(-x)
  }
}

///|
pub impl Density for Exponential with fn logpdf(_self, x) {
  if x < 0.0 {
    neg_inf()
  } else {
    -x
  }
}

///|
pub impl Cdf for Exponential with fn cdf(_self, x) {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-x)
  }
}

///|
pub impl Cdf for Exponential with fn sf(_self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-x)
  }
}

///|
pub impl Quantile for Exponential with fn quantile(_self, p) {
  -ln1p(-p)
}

///|
pub impl Quantile for Exponential with fn isf(_self, q) {
  -@math.ln(q)
}

///|
pub impl Density for Laplace with fn pdf(_self, x) {
  0.5 * @math.exp(-x.abs())
}

///|
pub impl Density for Laplace with fn logpdf(_self, x) {
  -@math.ln(2.0) - x.abs()
}

///|
pub impl Cdf for Laplace with fn cdf(_self, x) {
  if x < 0.0 {
    0.5 * @math.exp(x)
  } else {
    1.0 - 0.5 * @math.exp(-x)
  }
}

///|
pub impl Cdf for Laplace with fn sf(_self, x) {
  if x > 0.0 {
    0.5 * @math.exp(-x)
  } else {
    1.0 - 0.5 * @math.exp(x)
  }
}

///|
pub impl Quantile for Laplace with fn quantile(_self, p) {
  if p <= 0.5 {
    @math.ln(2.0 * p)
  } else {
    -@math.ln(2.0 * (1.0 - p))
  }
}

///|
pub impl Quantile for Laplace with fn isf(self, q) {
  -Quantile::quantile(self, q)
}

///|
pub impl Density for Logistic with fn pdf(self, x) {
  @math.exp(Density::logpdf(self, x))
}

///|
pub impl Density for Logistic with fn logpdf(_self, x) {
  let a = x.abs()
  -a - 2.0 * ln1p(@math.exp(-a))
}

///|
pub impl Cdf for Logistic with fn cdf(_self, x) {
  if x >= 0.0 {
    1.0 / (1.0 + @math.exp(-x))
  } else {
    let e = @math.exp(x)
    e / (1.0 + e)
  }
}

///|
pub impl Cdf for Logistic with fn sf(self, x) {
  Cdf::cdf(self, -x)
}

///|
pub impl Quantile for Logistic with fn quantile(_self, p) {
  @math.ln(p) - ln1p(-p)
}

///|
pub impl Quantile for Logistic with fn isf(self, q) {
  -Quantile::quantile(self, q)
}

///|
pub impl Density for Cauchy with fn pdf(_self, x) {
  1.0 / (@math.PI * (1.0 + x * x))
}

///|
/// 1/2 + atan(x)/π, written as atan(-1/x)/π for x < 0 so the lower tail
/// keeps its digits.
pub impl Cdf for Cauchy with fn cdf(_self, x) {
  if x < 0.0 {
    @math.atan(-1.0 / x) / @math.PI
  } else if x == 0.0 {
    0.5
  } else {
    1.0 - @math.atan(1.0 / x) / @math.PI
  }
}

///|
pub impl Cdf for Cauchy with fn sf(self, x) {
  Cdf::cdf(self, -x)
}

///|
pub impl Quantile for Cauchy with fn quantile(_self, p) {
  if p < 0.5 {
    -1.0 / @math.tan(@math.PI * p)
  } else if p == 0.5 {
    0.0
  } else {
    1.0 / @math.tan(@math.PI * (1.0 - p))
  }
}

///|
pub impl Quantile for Cauchy with fn isf(self, q) {
  -Quantile::quantile(self, q)
}

// ---------------------------------------------------------------------------
// Gumbel, Weibull, Pareto, Lomax, Rayleigh
// ---------------------------------------------------------------------------

///|
pub impl Density for GumbelRight with fn pdf(self, x) {
  @math.exp(Density::logpdf(self, x))
}

///|
pub impl Density for GumbelRight with fn logpdf(_self, x) {
  -x - @math.exp(-x)
}

///|
pub impl Cdf for GumbelRight with fn cdf(_self, x) {
  @math.exp(-@math.exp(-x))
}

///|
pub impl Cdf for GumbelRight with fn sf(_self, x) {
  -@math.expm1(-@math.exp(-x))
}

///|
pub impl Quantile for GumbelRight with fn quantile(_self, p) {
  -@math.ln(-@math.ln(p))
}

///|
pub impl Quantile for GumbelRight with fn isf(_self, q) {
  -@math.ln(-ln1p(-q))
}

///|
pub impl Density for GumbelLeft with fn pdf(self, x) {
  @math.exp(Density::logpdf(self, x))
}

///|
pub impl Density for GumbelLeft with fn logpdf(_self, x) {
  x - @math.exp(x)
}

///|
pub impl Cdf for GumbelLeft with fn cdf(_self, x) {
  -@math.expm1(-@math.exp(x))
}

///|
pub impl Cdf for GumbelLeft with fn sf(_self, x) {
  @math.exp(-@math.exp(x))
}

///|
pub impl Quantile for GumbelLeft with fn quantile(_self, p) {
  @math.ln(-ln1p(-p))
}

///|
pub impl Quantile for GumbelLeft with fn isf(_self, q) {
  @math.ln(-@math.ln(q))
}

///|
pub impl Density for Weibull with fn pdf(self, x) {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Weibull with fn logpdf(self, x) {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) + (self.c - 1.0) * @math.ln(x) - @math.pow(x, self.c)
  }
}

///|
pub impl Cdf for Weibull with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-@math.pow(x, self.c))
  }
}

///|
pub impl Cdf for Weibull with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-@math.pow(x, self.c))
  }
}

///|
pub impl Quantile for Weibull with fn quantile(self, p) {
  @math.pow(-ln1p(-p), 1.0 / self.c)
}

///|
pub impl Quantile for Weibull with fn isf(self, q) {
  @math.pow(-@math.ln(q), 1.0 / self.c)
}

///|
pub impl Density for WeibullMax with fn pdf(self, x) {
  if x > 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for WeibullMax with fn logpdf(self, x) {
  if x > 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) + (self.c - 1.0) * @math.ln(-x) - @math.pow(-x, self.c)
  }
}

///|
pub impl Cdf for WeibullMax with fn cdf(self, x) {
  if x >= 0.0 {
    1.0
  } else {
    @math.exp(-@math.pow(-x, self.c))
  }
}

///|
pub impl Cdf for WeibullMax with fn sf(self, x) {
  if x >= 0.0 {
    0.0
  } else {
    -@math.expm1(-@math.pow(-x, self.c))
  }
}

///|
pub impl Quantile for WeibullMax with fn quantile(self, p) {
  -@math.pow(-@math.ln(p), 1.0 / self.c)
}

///|
pub impl Quantile for WeibullMax with fn isf(self, q) {
  -@math.pow(-ln1p(-q), 1.0 / self.c)
}

///|
pub impl Density for Pareto with fn pdf(self, x) {
  if x < 1.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Pareto with fn logpdf(self, x) {
  if x < 1.0 {
    neg_inf()
  } else {
    @math.ln(self.b) - (self.b + 1.0) * @math.ln(x)
  }
}

///|
pub impl Cdf for Pareto with fn cdf(self, x) {
  if x <= 1.0 {
    0.0
  } else {
    -@math.expm1(-self.b * @math.ln(x))
  }
}

///|
pub impl Cdf for Pareto with fn sf(self, x) {
  if x <= 1.0 {
    1.0
  } else {
    @math.exp(-self.b * @math.ln(x))
  }
}

///|
pub impl Quantile for Pareto with fn quantile(self, p) {
  @math.exp(-ln1p(-p) / self.b)
}

///|
pub impl Quantile for Pareto with fn isf(self, q) {
  @math.exp(-@math.ln(q) / self.b)
}

///|
pub impl Density for Lomax with fn pdf(self, x) {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Lomax with fn logpdf(self, x) {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) - (self.c + 1.0) * ln1p(x)
  }
}

///|
pub impl Cdf for Lomax with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-self.c * ln1p(x))
  }
}

///|
pub impl Cdf for Lomax with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-self.c * ln1p(x))
  }
}

///|
pub impl Quantile for Lomax with fn quantile(self, p) {
  @math.expm1(-ln1p(-p) / self.c)
}

///|
pub impl Quantile for Lomax with fn isf(self, q) {
  @math.expm1(-@math.ln(q) / self.c)
}

///|
pub impl Density for Rayleigh with fn pdf(_self, x) {
  if x < 0.0 {
    0.0
  } else {
    x * @math.exp(-0.5 * x * x)
  }
}

///|
pub impl Density for Rayleigh with fn logpdf(_self, x) {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(x) - 0.5 * x * x
  }
}

///|
pub impl Cdf for Rayleigh with fn cdf(_self, x) {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-0.5 * x * x)
  }
}

///|
pub impl Cdf for Rayleigh with fn sf(_self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-0.5 * x * x)
  }
}

///|
pub impl Quantile for Rayleigh with fn quantile(_self, p) {
  (-2.0 * ln1p(-p)).sqrt()
}

///|
pub impl Quantile for Rayleigh with fn isf(_self, q) {
  (-2.0 * @math.ln(q)).sqrt()
}

// ---------------------------------------------------------------------------
// Bounded: triangular, arcsine, power law, Bradford, log-uniform, anglit,
// truncated exponential, truncated Pareto
// ---------------------------------------------------------------------------

///|
pub impl Density for Triangular with fn pdf(self, x) {
  if x < 0.0 || x > 1.0 {
    0.0
  } else if x < self.c {
    2.0 * x / self.c
  } else {
    2.0 * (1.0 - x) / (1.0 - self.c)
  }
}

///|
pub impl Cdf for Triangular with fn cdf(self, x) {
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
pub impl Cdf for Triangular with fn sf(self, x) {
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
pub impl Quantile for Triangular with fn quantile(self, p) {
  if p < self.c {
    (self.c * p).sqrt()
  } else {
    1.0 - ((1.0 - self.c) * (1.0 - p)).sqrt()
  }
}

///|
pub impl Quantile for Triangular with fn isf(self, q) {
  if q > 1.0 - self.c {
    (self.c * (1.0 - q)).sqrt()
  } else {
    1.0 - ((1.0 - self.c) * q).sqrt()
  }
}

///|
pub impl Density for Arcsine with fn pdf(_self, x) {
  if x <= 0.0 || x >= 1.0 {
    0.0
  } else {
    1.0 / (@math.PI * (x * (1.0 - x)).sqrt())
  }
}

///|
pub impl Cdf for Arcsine with fn cdf(_self, x) {
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
pub impl Cdf for Arcsine with fn sf(_self, x) {
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
pub impl Quantile for Arcsine with fn quantile(_self, p) {
  let s = @math.sin(0.5 * @math.PI * p)
  s * s
}

///|
pub impl Quantile for Arcsine with fn isf(_self, q) {
  let c = @math.cos(0.5 * @math.PI * q)
  c * c
}

///|
pub impl Density for PowerLaw with fn pdf(self, x) {
  if x < 0.0 || x > 1.0 {
    0.0
  } else {
    self.a * @math.pow(x, self.a - 1.0)
  }
}

///|
pub impl Cdf for PowerLaw with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    @math.pow(x, self.a)
  }
}

///|
pub impl Cdf for PowerLaw with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else {
    -@math.expm1(self.a * @math.ln(x))
  }
}

///|
pub impl Quantile for PowerLaw with fn quantile(self, p) {
  @math.pow(p, 1.0 / self.a)
}

///|
pub impl Quantile for PowerLaw with fn isf(self, q) {
  @math.exp(ln1p(-q) / self.a)
}

///|
pub impl Density for Bradford with fn pdf(self, x) {
  if x < 0.0 || x > 1.0 {
    0.0
  } else {
    self.c / (ln1p(self.c) * (1.0 + self.c * x))
  }
}

///|
pub impl Cdf for Bradford with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    ln1p(self.c * x) / ln1p(self.c)
  }
}

///|
pub impl Cdf for Bradford with fn sf(self, x) {
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
pub impl Quantile for Bradford with fn quantile(self, p) {
  @math.expm1(p * ln1p(self.c)) / self.c
}

///|
pub impl Quantile for Bradford with fn isf(self, q) {
  @math.expm1((1.0 - q) * ln1p(self.c)) / self.c
}

///|
pub impl Density for LogUniform with fn pdf(self, x) {
  if x < self.a || x > self.b {
    0.0
  } else {
    1.0 / (x * (@math.ln(self.b) - @math.ln(self.a)))
  }
}

///|
pub impl Cdf for LogUniform with fn cdf(self, x) {
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
pub impl Cdf for LogUniform with fn sf(self, x) {
  if x <= self.a {
    1.0
  } else if x >= self.b {
    0.0
  } else {
    ln1p((self.b - x) / x) / (@math.ln(self.b) - @math.ln(self.a))
  }
}

///|
pub impl Quantile for LogUniform with fn quantile(self, p) {
  @math.exp(@math.ln(self.a) + p * (@math.ln(self.b) - @math.ln(self.a)))
}

///|
pub impl Quantile for LogUniform with fn isf(self, q) {
  @math.exp(@math.ln(self.b) - q * (@math.ln(self.b) - @math.ln(self.a)))
}

///|
pub impl Density for Anglit with fn pdf(_self, x) {
  if x < -0.25 * @math.PI || x > 0.25 * @math.PI {
    0.0
  } else {
    @math.cos(2.0 * x)
  }
}

///|
pub impl Cdf for Anglit with fn cdf(_self, x) {
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
pub impl Cdf for Anglit with fn sf(_self, x) {
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
pub impl Quantile for Anglit with fn quantile(_self, p) {
  @math.asin(p.sqrt()) - 0.25 * @math.PI
}

///|
pub impl Quantile for Anglit with fn isf(_self, q) {
  @math.acos(q.sqrt()) - 0.25 * @math.PI
}

///|
pub impl Density for TruncatedExponential with fn pdf(self, x) {
  if x < 0.0 || x > self.b {
    0.0
  } else {
    @math.exp(-x) / -@math.expm1(-self.b)
  }
}

///|
pub impl Cdf for TruncatedExponential with fn cdf(self, x) {
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
pub impl Cdf for TruncatedExponential with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else if x >= self.b {
    0.0
  } else {
    @math.exp(-x) * @math.expm1(x - self.b) / @math.expm1(-self.b)
  }
}

///|
pub impl Quantile for TruncatedExponential with fn quantile(self, p) {
  -ln1p(p * @math.expm1(-self.b))
}

///|
pub impl Quantile for TruncatedExponential with fn isf(self, q) {
  -@math.ln(@math.exp(-self.b) - q * @math.expm1(-self.b))
}

///|
fn TruncatedPareto::mass(self : TruncatedPareto) -> Double {
  -@math.expm1(-self.b * @math.ln(self.c))
}

///|
pub impl Density for TruncatedPareto with fn pdf(self, x) {
  if x < 1.0 || x > self.c {
    0.0
  } else {
    self.b * @math.exp(-(self.b + 1.0) * @math.ln(x)) / self.mass()
  }
}

///|
pub impl Cdf for TruncatedPareto with fn cdf(self, x) {
  if x <= 1.0 {
    0.0
  } else if x >= self.c {
    1.0
  } else {
    -@math.expm1(-self.b * @math.ln(x)) / self.mass()
  }
}

///|
pub impl Cdf for TruncatedPareto with fn sf(self, x) {
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
pub impl Quantile for TruncatedPareto with fn quantile(self, p) {
  @math.exp(-ln1p(-p * self.mass()) / self.b)
}

///|
pub impl Quantile for TruncatedPareto with fn isf(self, q) {
  @math.exp(
    -@math.ln(@math.exp(-self.b * @math.ln(self.c)) + q * self.mass()) / self.b,
  )
}

// ---------------------------------------------------------------------------
// Half and folded forms, hyperbolic secant
// ---------------------------------------------------------------------------

///|
pub impl Density for HalfCauchy with fn pdf(_self, x) {
  if x < 0.0 {
    0.0
  } else {
    2.0 / (@math.PI * (1.0 + x * x))
  }
}

///|
pub impl Cdf for HalfCauchy with fn cdf(_self, x) {
  if x <= 0.0 {
    0.0
  } else {
    2.0 / @math.PI * @math.atan(x)
  }
}

///|
pub impl Cdf for HalfCauchy with fn sf(_self, x) {
  if x <= 0.0 {
    1.0
  } else {
    2.0 / @math.PI * @math.atan(1.0 / x)
  }
}

///|
/// tan(πp/2) is ill-conditioned near p = 1; use 1/tan(π(1 - p)/2) there.
pub impl Quantile for HalfCauchy with fn quantile(_self, p) {
  if p <= 0.5 {
    @math.tan(0.5 * @math.PI * p)
  } else {
    1.0 / @math.tan(0.5 * @math.PI * (1.0 - p))
  }
}

///|
pub impl Quantile for HalfCauchy with fn isf(_self, q) {
  1.0 / @math.tan(0.5 * @math.PI * q)
}

///|
pub impl Density for HalfLogistic with fn pdf(self, x) {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for HalfLogistic with fn logpdf(_self, x) {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(2.0) - x - 2.0 * ln1p(@math.exp(-x))
  }
}

///|
pub impl Cdf for HalfLogistic with fn cdf(_self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.tanh(0.5 * x)
  }
}

///|
pub impl Cdf for HalfLogistic with fn sf(_self, x) {
  if x <= 0.0 {
    1.0
  } else {
    let e = @math.exp(-x)
    2.0 * e / (1.0 + e)
  }
}

///|
pub impl Quantile for HalfLogistic with fn quantile(_self, p) {
  ln1p(p) - ln1p(-p)
}

///|
/// ln((2 - q)/q) = ln1p(2(1 - q)/q), exact when q is near 1.
pub impl Quantile for HalfLogistic with fn isf(_self, q) {
  ln1p(2.0 * (1.0 - q) / q)
}

///|
pub impl Density for HyperbolicSecant with fn pdf(_self, x) {
  1.0 / (@math.PI * @math.cosh(x))
}

///|
pub impl Cdf for HyperbolicSecant with fn cdf(_self, x) {
  2.0 / @math.PI * @math.atan(@math.exp(x))
}

///|
pub impl Cdf for HyperbolicSecant with fn sf(_self, x) {
  2.0 / @math.PI * @math.atan(@math.exp(-x))
}

///|
pub impl Quantile for HyperbolicSecant with fn quantile(_self, p) {
  if p <= 0.5 {
    @math.ln(@math.tan(0.5 * @math.PI * p))
  } else {
    -@math.ln(@math.tan(0.5 * @math.PI * (1.0 - p)))
  }
}

///|
pub impl Quantile for HyperbolicSecant with fn isf(self, q) {
  -Quantile::quantile(self, q)
}

// ---------------------------------------------------------------------------
// Generalised logistic, Pareto, extreme value
// ---------------------------------------------------------------------------

///|
pub impl Density for GeneralizedLogistic with fn pdf(self, x) {
  @math.exp(Density::logpdf(self, x))
}

///|
pub impl Density for GeneralizedLogistic with fn logpdf(self, x) {
  @math.ln(self.c) - x - (self.c + 1.0) * @num.log1p_exp(-x)
}

///|
pub impl Cdf for GeneralizedLogistic with fn cdf(self, x) {
  @math.exp(-self.c * @num.log1p_exp(-x))
}

///|
pub impl Cdf for GeneralizedLogistic with fn sf(self, x) {
  -@math.expm1(-self.c * @num.log1p_exp(-x))
}

///|
pub impl Quantile for GeneralizedLogistic with fn quantile(self, p) {
  -@math.ln(@math.expm1(-@math.ln(p) / self.c))
}

///|
pub impl Quantile for GeneralizedLogistic with fn isf(self, q) {
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
pub impl Density for GeneralizedPareto with fn pdf(self, x) {
  if x < 0.0 || x > self.upper() {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for GeneralizedPareto with fn logpdf(self, x) {
  if x < 0.0 || x > self.upper() {
    neg_inf()
  } else if self.c == 0.0 {
    -x
  } else {
    -(1.0 + 1.0 / self.c) * ln1p(self.c * x)
  }
}

///|
pub impl Cdf for GeneralizedPareto with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else if x >= self.upper() {
    1.0
  } else {
    -@math.expm1(self.log_sf(x))
  }
}

///|
pub impl Cdf for GeneralizedPareto with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else if x >= self.upper() {
    0.0
  } else {
    @math.exp(self.log_sf(x))
  }
}

///|
pub impl Quantile for GeneralizedPareto with fn quantile(self, p) {
  if self.c == 0.0 {
    -ln1p(-p)
  } else {
    @math.expm1(-self.c * ln1p(-p)) / self.c
  }
}

///|
pub impl Quantile for GeneralizedPareto with fn isf(self, q) {
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
pub impl Density for GeneralizedExtremeValue with fn pdf(self, x) {
  if !self.inside(x) {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for GeneralizedExtremeValue with fn logpdf(self, x) {
  if !self.inside(x) {
    neg_inf()
  } else if self.c == 0.0 {
    -x - @math.exp(-x)
  } else {
    -self.t(x) + (1.0 / self.c - 1.0) * ln1p(-self.c * x)
  }
}

///|
pub impl Cdf for GeneralizedExtremeValue with fn cdf(self, x) {
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
pub impl Cdf for GeneralizedExtremeValue with fn sf(self, x) {
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
pub impl Quantile for GeneralizedExtremeValue with fn quantile(self, p) {
  let l = -@math.ln(p)
  if self.c == 0.0 {
    -@math.ln(l)
  } else {
    -@math.expm1(self.c * @math.ln(l)) / self.c
  }
}

///|
pub impl Quantile for GeneralizedExtremeValue with fn isf(self, q) {
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
pub impl Density for Burr with fn pdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Burr with fn logpdf(self, x) {
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
pub impl Cdf for Burr with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(burr3_log_cdf(self.c, self.d, x))
  }
}

///|
pub impl Cdf for Burr with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    -@math.expm1(burr3_log_cdf(self.c, self.d, x))
  }
}

///|
pub impl Quantile for Burr with fn quantile(self, p) {
  @math.pow(@math.expm1(-@math.ln(p) / self.d), -1.0 / self.c)
}

///|
pub impl Quantile for Burr with fn isf(self, q) {
  @math.pow(@math.expm1(-ln1p(-q) / self.d), -1.0 / self.c)
}

///|
pub impl Density for Fisk with fn pdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Fisk with fn logpdf(self, x) {
  if x <= 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) -
    (self.c + 1.0) * @math.ln(x) -
    2.0 * ln1p(@math.pow(x, -self.c))
  }
}

///|
pub impl Cdf for Fisk with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(burr3_log_cdf(self.c, 1.0, x))
  }
}

///|
pub impl Cdf for Fisk with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    -@math.expm1(burr3_log_cdf(self.c, 1.0, x))
  }
}

///|
pub impl Quantile for Fisk with fn quantile(self, p) {
  @math.pow(@math.expm1(-@math.ln(p)), -1.0 / self.c)
}

///|
pub impl Quantile for Fisk with fn isf(self, q) {
  @math.pow(@math.expm1(-ln1p(-q)), -1.0 / self.c)
}

///|
pub impl Density for Burr12 with fn pdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Burr12 with fn logpdf(self, x) {
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
pub impl Cdf for Burr12 with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-self.d * ln1p(@math.pow(x, self.c)))
  }
}

///|
pub impl Cdf for Burr12 with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-self.d * ln1p(@math.pow(x, self.c)))
  }
}

///|
pub impl Quantile for Burr12 with fn quantile(self, p) {
  @math.pow(@math.expm1(-ln1p(-p) / self.d), 1.0 / self.c)
}

///|
pub impl Quantile for Burr12 with fn isf(self, q) {
  @math.pow(@math.expm1(-@math.ln(q) / self.d), 1.0 / self.c)
}

///|
/// F = x (a + x^a)^(-1/a) = exp(-ln(1 + a x^-a)/a).
pub impl Density for Kappa3 with fn pdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    self.a * @math.pow(self.a + @math.pow(x, self.a), -1.0 / self.a - 1.0)
  }
}

///|
pub impl Cdf for Kappa3 with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(-ln1p(self.a * @math.pow(x, -self.a)) / self.a)
  }
}

///|
pub impl Cdf for Kappa3 with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    -@math.expm1(-ln1p(self.a * @math.pow(x, -self.a)) / self.a)
  }
}

///|
pub impl Quantile for Kappa3 with fn quantile(self, p) {
  @math.pow(self.a / @math.expm1(-self.a * @math.ln(p)), 1.0 / self.a)
}

///|
pub impl Quantile for Kappa3 with fn isf(self, q) {
  @math.pow(self.a / @math.expm1(-self.a * ln1p(-q)), 1.0 / self.a)
}

// ---------------------------------------------------------------------------
// Exponential variants, Gompertz, double Weibull, log-Laplace,
// asymmetric Laplace, skew Cauchy
// ---------------------------------------------------------------------------

///|
pub impl Density for ExponentialPower with fn pdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for ExponentialPower with fn logpdf(self, x) {
  if x <= 0.0 {
    neg_inf()
  } else {
    let xb = @math.pow(x, self.b)
    @math.ln(self.b) + (self.b - 1.0) * @math.ln(x) + xb - @math.expm1(xb)
  }
}

///|
pub impl Cdf for ExponentialPower with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-@math.expm1(@math.pow(x, self.b)))
  }
}

///|
pub impl Cdf for ExponentialPower with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-@math.expm1(@math.pow(x, self.b)))
  }
}

///|
pub impl Quantile for ExponentialPower with fn quantile(self, p) {
  @math.pow(ln1p(-ln1p(-p)), 1.0 / self.b)
}

///|
pub impl Quantile for ExponentialPower with fn isf(self, q) {
  @math.pow(ln1p(-@math.ln(q)), 1.0 / self.b)
}

///|
pub impl Density for ExponentiatedWeibull with fn pdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for ExponentiatedWeibull with fn logpdf(self, x) {
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
pub impl Cdf for ExponentiatedWeibull with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @math.exp(self.a * @math.ln(-@math.expm1(-@math.pow(x, self.c))))
  }
}

///|
pub impl Cdf for ExponentiatedWeibull with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    -@math.expm1(self.a * @math.ln(-@math.expm1(-@math.pow(x, self.c))))
  }
}

///|
/// x^c = -ln(1 - p^(1/a)), with 1 - p^(1/a) = -expm1(ln(p)/a) so that
/// p^(1/a) near 1 does not cancel.
pub impl Quantile for ExponentiatedWeibull with fn quantile(self, p) {
  @math.pow(-@math.ln(-@math.expm1(@math.ln(p) / self.a)), 1.0 / self.c)
}

///|
pub impl Quantile for ExponentiatedWeibull with fn isf(self, q) {
  @math.pow(-@math.ln(-@math.expm1(ln1p(-q) / self.a)), 1.0 / self.c)
}

///|
pub impl Density for Gompertz with fn pdf(self, x) {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Gompertz with fn logpdf(self, x) {
  if x < 0.0 {
    neg_inf()
  } else {
    @math.ln(self.c) + x - self.c * @math.expm1(x)
  }
}

///|
pub impl Cdf for Gompertz with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    -@math.expm1(-self.c * @math.expm1(x))
  }
}

///|
pub impl Cdf for Gompertz with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @math.exp(-self.c * @math.expm1(x))
  }
}

///|
pub impl Quantile for Gompertz with fn quantile(self, p) {
  ln1p(-ln1p(-p) / self.c)
}

///|
pub impl Quantile for Gompertz with fn isf(self, q) {
  ln1p(-@math.ln(q) / self.c)
}

///|
pub impl Density for DoubleWeibull with fn pdf(self, x) {
  let a = x.abs()
  0.5 * self.c * @math.pow(a, self.c - 1.0) * @math.exp(-@math.pow(a, self.c))
}

///|
pub impl Cdf for DoubleWeibull with fn cdf(self, x) {
  if x <= 0.0 {
    0.5 * @math.exp(-@math.pow(-x, self.c))
  } else {
    1.0 - 0.5 * @math.exp(-@math.pow(x, self.c))
  }
}

///|
pub impl Cdf for DoubleWeibull with fn sf(self, x) {
  Cdf::cdf(self, -x)
}

///|
pub impl Quantile for DoubleWeibull with fn quantile(self, p) {
  if p <= 0.5 {
    -@math.pow(-@math.ln(2.0 * p), 1.0 / self.c)
  } else {
    @math.pow(-@math.ln(2.0 * (1.0 - p)), 1.0 / self.c)
  }
}

///|
pub impl Quantile for DoubleWeibull with fn isf(self, q) {
  -Quantile::quantile(self, q)
}

///|
pub impl Density for LogLaplace with fn pdf(self, x) {
  if x <= 0.0 {
    0.0
  } else if x < 1.0 {
    0.5 * self.c * @math.pow(x, self.c - 1.0)
  } else {
    0.5 * self.c * @math.pow(x, -self.c - 1.0)
  }
}

///|
pub impl Cdf for LogLaplace with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else if x < 1.0 {
    0.5 * @math.pow(x, self.c)
  } else {
    1.0 - 0.5 * @math.pow(x, -self.c)
  }
}

///|
pub impl Cdf for LogLaplace with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else if x < 1.0 {
    1.0 - 0.5 * @math.pow(x, self.c)
  } else {
    0.5 * @math.pow(x, -self.c)
  }
}

///|
pub impl Quantile for LogLaplace with fn quantile(self, p) {
  if p < 0.5 {
    @math.pow(2.0 * p, 1.0 / self.c)
  } else {
    @math.pow(2.0 * (1.0 - p), -1.0 / self.c)
  }
}

///|
pub impl Quantile for LogLaplace with fn isf(self, q) {
  if q > 0.5 {
    @math.pow(2.0 * (1.0 - q), 1.0 / self.c)
  } else {
    @math.pow(2.0 * q, -1.0 / self.c)
  }
}

///|
pub impl Density for AsymmetricLaplace with fn pdf(self, x) {
  let k = self.kappa
  let e = if x >= 0.0 { @math.exp(-x * k) } else { @math.exp(x / k) }
  e / (k + 1.0 / k)
}

///|
pub impl Cdf for AsymmetricLaplace with fn cdf(self, x) {
  let k2 = self.kappa * self.kappa
  if x >= 0.0 {
    1.0 - @math.exp(-x * self.kappa) / (1.0 + k2)
  } else {
    k2 * @math.exp(x / self.kappa) / (1.0 + k2)
  }
}

///|
pub impl Cdf for AsymmetricLaplace with fn sf(self, x) {
  let k2 = self.kappa * self.kappa
  if x >= 0.0 {
    @math.exp(-x * self.kappa) / (1.0 + k2)
  } else {
    1.0 - k2 * @math.exp(x / self.kappa) / (1.0 + k2)
  }
}

///|
pub impl Quantile for AsymmetricLaplace with fn quantile(self, p) {
  let k = self.kappa
  let k2 = k * k
  if p <= k2 / (1.0 + k2) {
    k * @math.ln(p * (1.0 + k2) / k2)
  } else {
    -@math.ln((1.0 - p) * (1.0 + k2)) / k
  }
}

///|
pub impl Quantile for AsymmetricLaplace with fn isf(self, q) {
  let k = self.kappa
  let k2 = k * k
  if q <= 1.0 / (1.0 + k2) {
    -@math.ln(q * (1.0 + k2)) / k
  } else {
    k * @math.ln((1.0 - q) * (1.0 + k2) / k2)
  }
}

///|
pub impl Density for SkewCauchy with fn pdf(self, x) {
  let s = if x < 0.0 { 1.0 - self.a } else { 1.0 + self.a }
  1.0 / (@math.PI * (x * x / (s * s) + 1.0))
}

///|
pub impl Cdf for SkewCauchy with fn cdf(self, x) {
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
pub impl Cdf for SkewCauchy with fn sf(self, x) {
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
pub impl Quantile for SkewCauchy with fn quantile(self, p) {
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
pub impl Quantile for SkewCauchy with fn isf(self, q) {
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
// Discrete
// ---------------------------------------------------------------------------

///|
/// The smallest k in [lo, hi] with P(X <= k) >= q, starting from a guess
/// from a continuous inversion and correcting it by exact comparisons.
fn[D : IntCdf] refine_quantile(
  d : D,
  q : Double,
  guess : Int32,
  lo : Int32,
  hi : Int32,
) -> Int32 {
  let mut k = if guess < lo { lo } else if guess > hi { hi } else { guess }
  while k > lo && IntCdf::cdf_int(d, k - 1) >= q {
    k = k - 1
  }
  while k < hi && IntCdf::cdf_int(d, k) < q {
    k = k + 1
  }
  k
}

///|
fn ceil_int(x : Double) -> Int32 {
  if x.is_nan() {
    0
  } else if x > 2.0e9 {
    2147483647
  } else if x < -2.0e9 {
    -2147483647
  } else {
    @math.ceil(x).to_int()
  }
}

///|
pub impl Mass for Bernoulli with fn pmf(self, k) {
  if k == 0 {
    1.0 - self.p
  } else if k == 1 {
    self.p
  } else {
    0.0
  }
}

///|
pub impl IntCdf for Bernoulli with fn cdf_int(self, k) {
  if k < 0 {
    0.0
  } else if k < 1 {
    1.0 - self.p
  } else {
    1.0
  }
}

///|
pub impl IntCdf for Bernoulli with fn sf_int(self, k) {
  if k < 0 {
    1.0
  } else if k < 1 {
    self.p
  } else {
    0.0
  }
}

///|
pub impl IntQuantile for Bernoulli with fn quantile_int(self, q) {
  if q <= 1.0 - self.p {
    0
  } else {
    1
  }
}

///|
pub impl Mass for Geometric with fn pmf(self, k) {
  if k < 1 {
    0.0
  } else {
    self.p * @math.exp((k - 1).to_double() * ln1p(-self.p))
  }
}

///|
pub impl IntCdf for Geometric with fn cdf_int(self, k) {
  if k < 1 {
    0.0
  } else {
    -@math.expm1(k.to_double() * ln1p(-self.p))
  }
}

///|
pub impl IntCdf for Geometric with fn sf_int(self, k) {
  if k < 1 {
    1.0
  } else {
    @math.exp(k.to_double() * ln1p(-self.p))
  }
}

///|
pub impl IntQuantile for Geometric with fn quantile_int(self, q) {
  let guess = ceil_int(ln1p(-q) / ln1p(-self.p))
  refine_quantile(self, q, guess, 1, 2147483647)
}

///|
pub impl Mass for DiscreteUniform with fn pmf(self, k) {
  if k < self.low || k >= self.high {
    0.0
  } else {
    1.0 / (self.high - self.low).to_double()
  }
}

///|
pub impl IntCdf for DiscreteUniform with fn cdf_int(self, k) {
  if k < self.low {
    0.0
  } else if k >= self.high - 1 {
    1.0
  } else {
    (k - self.low + 1).to_double() / (self.high - self.low).to_double()
  }
}

///|
pub impl IntCdf for DiscreteUniform with fn sf_int(self, k) {
  if k < self.low {
    1.0
  } else if k >= self.high - 1 {
    0.0
  } else {
    (self.high - 1 - k).to_double() / (self.high - self.low).to_double()
  }
}

///|
pub impl IntQuantile for DiscreteUniform with fn quantile_int(self, q) {
  let n = (self.high - self.low).to_double()
  refine_quantile(
    self,
    q,
    ceil_int(q * n) + self.low - 1,
    self.low,
    self.high - 1,
  )
}

///|
pub impl Mass for DiscreteLaplace with fn pmf(self, k) {
  @math.tanh(0.5 * self.a) * @math.exp(-self.a * k.abs().to_double())
}

///|
pub impl IntCdf for DiscreteLaplace with fn cdf_int(self, k) {
  let z = 1.0 + @math.exp(-self.a)
  if k >= 0 {
    1.0 - @math.exp(-self.a * (k + 1).to_double()) / z
  } else {
    @math.exp(self.a * k.to_double()) / z
  }
}

///|
pub impl IntCdf for DiscreteLaplace with fn sf_int(self, k) {
  let z = 1.0 + @math.exp(-self.a)
  if k >= 0 {
    @math.exp(-self.a * (k + 1).to_double()) / z
  } else {
    1.0 - @math.exp(self.a * k.to_double()) / z
  }
}

///|
pub impl IntQuantile for DiscreteLaplace with fn quantile_int(self, q) {
  let z = 1.0 + @math.exp(-self.a)
  let guess = if q < 1.0 / z {
    @math.floor(@math.ln(q * z) / self.a).to_int()
  } else {
    ceil_int(-@math.ln((1.0 - q) * z) / self.a - 1.0)
  }
  refine_quantile(self, q, guess, -2147483647, 2147483647)
}

///|
pub impl Mass for Planck with fn pmf(self, k) {
  if k < 0 {
    0.0
  } else {
    -@math.expm1(-self.lambda) * @math.exp(-self.lambda * k.to_double())
  }
}

///|
pub impl IntCdf for Planck with fn cdf_int(self, k) {
  if k < 0 {
    0.0
  } else {
    -@math.expm1(-self.lambda * (k + 1).to_double())
  }
}

///|
pub impl IntCdf for Planck with fn sf_int(self, k) {
  if k < 0 {
    1.0
  } else {
    @math.exp(-self.lambda * (k + 1).to_double())
  }
}

///|
pub impl IntQuantile for Planck with fn quantile_int(self, q) {
  refine_quantile(
    self,
    q,
    ceil_int(-ln1p(-q) / self.lambda - 1.0),
    0,
    2147483647,
  )
}

///|
fn Boltzmann::z(self : Boltzmann) -> Double {
  -@math.expm1(-self.lambda * self.n.to_double())
}

///|
pub impl Mass for Boltzmann with fn pmf(self, k) {
  if k < 0 || k >= self.n {
    0.0
  } else {
    -@math.expm1(-self.lambda) *
    @math.exp(-self.lambda * k.to_double()) /
    self.z()
  }
}

///|
pub impl IntCdf for Boltzmann with fn cdf_int(self, k) {
  if k < 0 {
    0.0
  } else if k >= self.n - 1 {
    1.0
  } else {
    -@math.expm1(-self.lambda * (k + 1).to_double()) / self.z()
  }
}

///|
pub impl IntCdf for Boltzmann with fn sf_int(self, k) {
  if k < 0 {
    1.0
  } else if k >= self.n - 1 {
    0.0
  } else {
    @math.exp(-self.lambda * (k + 1).to_double()) *
    -@math.expm1(-self.lambda * (self.n - k - 1).to_double()) /
    self.z()
  }
}

///|
pub impl IntQuantile for Boltzmann with fn quantile_int(self, q) {
  let guess = ceil_int(-ln1p(-q * self.z()) / self.lambda - 1.0)
  refine_quantile(self, q, guess, 0, self.n - 1)
}

///|
fn Zipfian::norm(self : Zipfian) -> Double {
  let mut s = @num.CompensatedSum::zero()
  for j in 1..<=self.n {
    s = s.add(@math.pow(j.to_double(), -self.a))
  }
  s.value()
}

///|
pub impl Mass for Zipfian with fn pmf(self, k) {
  if k < 1 || k > self.n {
    0.0
  } else {
    @math.pow(k.to_double(), -self.a) / self.norm()
  }
}

///|
pub impl IntCdf for Zipfian with fn cdf_int(self, k) {
  if k < 1 {
    0.0
  } else if k >= self.n {
    1.0
  } else {
    let mut s = @num.CompensatedSum::zero()
    for j in 1..<=k {
      s = s.add(@math.pow(j.to_double(), -self.a))
    }
    s.value() / self.norm()
  }
}

///|
pub impl IntCdf for Zipfian with fn sf_int(self, k) {
  if k < 1 {
    1.0
  } else if k >= self.n {
    0.0
  } else {
    let mut s = @num.CompensatedSum::zero()
    for j in (k + 1)..<=self.n {
      s = s.add(@math.pow(j.to_double(), -self.a))
    }
    s.value() / self.norm()
  }
}

///|
pub impl IntQuantile for Zipfian with fn quantile_int(self, q) {
  refine_quantile(self, q, 1, 1, self.n)
}
