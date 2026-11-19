///|
/// Capabilities written by hand on top of the generated catalogue types:
/// the special-function paths of the reference distributions, moment and
/// characteristic functions, fitting, and alternative constructors. Field
/// names follow SciPy (generated in catalog.mbt).

// ---------------------------------------------------------------------------
// Normal: loc and scale built in.
// ---------------------------------------------------------------------------

///|
pub impl Density for Normal with fn pdf(self, x) {
  @math.exp(Density::logpdf(self, x))
}

///|
/// Closed form, so the far tails stay finite on the log scale.
pub impl Density for Normal with fn logpdf(self, x) {
  let z = (x - self.loc) / self.scale
  -0.5 * z * z - @math.ln(self.scale) - 0.5 * @math.ln(2.0 * @math.PI)
}

///|
pub impl Cdf for Normal with fn cdf(self, x) {
  @special.ndtr((x - self.loc) / self.scale)
}

///|
/// Symmetry gives the upper tail without cancellation.
pub impl Cdf for Normal with fn sf(self, x) {
  @special.ndtr((self.loc - x) / self.scale)
}

///|
pub impl Quantile for Normal with fn quantile(self, p) {
  self.loc + self.scale * @special.ndtri(p)
}

///|
/// By symmetry, never through 1 - q.
pub impl Quantile for Normal with fn isf(self, q) {
  self.loc - self.scale * @special.ndtri(q)
}

///|
pub impl Moments for Normal with fn mean(self) {
  Some(self.loc)
}

///|
pub impl Moments for Normal with fn variance(self) {
  Some(self.scale * self.scale)
}

///|
/// φ(t) = exp(i loc t - scale² t²/2).
pub impl Characteristic for Normal with fn cf(self, t) {
  let r = @math.exp(-0.5 * self.scale * self.scale * t * t)
  (r * @math.cos(self.loc * t), r * @math.sin(self.loc * t))
}

///|
/// M(t) = exp(loc t + scale² t²/2), finite for every t.
pub impl Mgf for Normal with fn mgf(self, t) {
  Some(@math.exp(self.loc * t + 0.5 * self.scale * self.scale * t * t))
}

///|
/// Maximum likelihood: the sample mean and the divisor-n standard
/// deviation. Needs two points, since one gives scale = 0.
pub impl Fit for Normal with fn fit(data) {
  if data.length() < 2 {
    raise TooFewPoints(needed=2, got=data.length())
  }
  ...
}

// ---------------------------------------------------------------------------
// Student t
// ---------------------------------------------------------------------------

///|
pub impl Moments for StudentT with fn mean(self) {
  if self.df > 1.0 {
    Some(0.0)
  } else {
    None
  }
}

///|
pub impl Moments for StudentT with fn variance(self) {
  if self.df > 2.0 {
    Some(self.df / (self.df - 2.0))
  } else if self.df > 1.0 {
    Some(1.0 / 0.0)
  } else {
    None
  }
}

///|
pub impl Density for StudentT with fn pdf(self, x) {
  @math.exp(Density::logpdf(self, x))
}

///|
/// -ln √ν - ln B(1/2, ν/2) - (ν+1)/2 · ln(1 + x²/ν).
pub impl Density for StudentT with fn logpdf(self, x) {
  -0.5 * @math.ln(self.df) -
  @special.lbeta(0.5, 0.5 * self.df) -
  0.5 * (self.df + 1.0) * @math.ln_1p(x * x / self.df)
}

///|
/// Half the lower tail, I_t(ν/2, 1/2)/2 with t = ν/(ν + x²); both t and
/// 1 - t = x²/(ν + x²) are formed directly, so neither end cancels.
fn StudentT::half_tail(self : StudentT, x : Double) -> Double {
  let d = self.df + x * x
  let (i, _) = @special.betainc_pair(0.5 * self.df, 0.5, self.df / d, x * x / d)
  0.5 * i
}

///|
pub impl Cdf for StudentT with fn cdf(self, x) {
  if x < 0.0 {
    self.half_tail(x)
  } else {
    1.0 - self.half_tail(x)
  }
}

///|
/// The p-value path of every t test, from the incomplete beta directly.
pub impl Cdf for StudentT with fn sf(self, x) {
  if x > 0.0 {
    self.half_tail(x)
  } else {
    1.0 - self.half_tail(x)
  }
}

///|
/// Inverts the tail it is in: I_t(ν/2, 1/2) = 2q for the far tail, and the
/// complementary I_z(1/2, ν/2) = 1 - 2q near the centre, so x stays
/// accurate both at 1e-300 and next to the median.
pub impl Quantile for StudentT with fn quantile(self, p) {
  if p == 0.5 {
    return 0.0
  }
  let q = if p < 0.5 { p } else { 1.0 - p }
  let target = 2.0 * q
  let x = if target <= 0.5 {
    let t = @special.betaincinv(0.5 * self.df, 0.5, target)
    (self.df * (1.0 - t) / t).sqrt()
  } else {
    let z = @special.betaincinv(0.5, 0.5 * self.df, 1.0 - target)
    (self.df * z / (1.0 - z)).sqrt()
  }
  if p < 0.5 {
    -x
  } else {
    x
  }
}

///|
pub impl Quantile for StudentT with fn isf(self, q) {
  -Quantile::quantile(self, q)
}

///|
/// No `Mgf`: E[exp(tX)] diverges for every t != 0. The characteristic
/// function is a modified Bessel function of the second kind.
pub impl Characteristic for StudentT with fn cf(self, t) {
  ignore((self, t))
  ...
}

// ---------------------------------------------------------------------------
// Gamma: shape a and scale built in (SciPy's parameterisation); the rate
// form is an alternative constructor.
// ---------------------------------------------------------------------------

///|
/// The same distribution by shape and rate (scale = 1/rate).
pub fn Gamma::with_rate(a : Double, rate : Double) -> Gamma raise ParamError {
  if !(rate > 0.0) || rate.is_inf() {
    raise InvalidParam(name="rate", value=rate, reason="must be in (0, inf)")
  }
  Gamma::new(a, scale=1.0 / rate)
}

///|
/// Gamma by shape and rate: a separate family value, so a fit reports its
/// parameters in the convention asked for.
pub fn Gamma::family_rate() -> Family[Gamma] {
  {
    name: "gamma (a, rate)",
    specs: [
      { name: "a", constraint: Positive, },
      { name: "rate", constraint: Positive, },
    ],
    make: p => Gamma::with_rate(p[0], p[1]),
    params: d => [d.a, 1.0 / d.scale],
  }
}

///|
pub impl Cdf for Gamma with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @special.gammainc(self.a, x / self.scale)
  }
}

///|
pub impl Cdf for Gamma with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @special.gammaincc(self.a, x / self.scale)
  }
}

///|
pub impl Density for Gamma with fn pdf(self, x) {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Gamma with fn logpdf(self, x) {
  if x < 0.0 {
    -1.0 / 0.0
  } else {
    let z = x / self.scale
    (self.a - 1.0) * @math.ln(z) -
    z -
    @special.lgamma(self.a) -
    @math.ln(self.scale)
  }
}

///|
pub impl Quantile for Gamma with fn quantile(self, p) {
  self.scale * @special.gammaincinv(self.a, p)
}

///|
pub impl Quantile for Gamma with fn isf(self, q) {
  self.scale * @special.gammainccinv(self.a, q)
}

///|
/// M(t) = (1 - t scale)^(-a) for t < 1/scale; it diverges from there on.
pub impl Mgf for Gamma with fn mgf(self, t) {
  if t * self.scale < 1.0 {
    Some(@math.pow(1.0 - t * self.scale, -self.a))
  } else {
    None
  }
}

// ---------------------------------------------------------------------------
// Binomial
// ---------------------------------------------------------------------------

///|
pub impl Mass for Binomial with fn pmf(self, k) {
  if k < 0 || k > self.n {
    0.0
  } else if self.p == 0.0 {
    if k == 0 {
      1.0
    } else {
      0.0
    }
  } else if self.p == 1.0 {
    if k == self.n {
      1.0
    } else {
      0.0
    }
  } else {
    @math.exp(
      @special.lchoose(self.n.to_double(), k.to_double()) +
      k.to_double() * @math.ln(self.p) +
      (self.n - k).to_double() * @math.ln_1p(-self.p),
    )
  }
}

///|
/// P(X <= k) = I_{1-p}(n - k, k + 1).
pub impl IntCdf for Binomial with fn cdf_int(self, k) {
  if k < 0 {
    0.0
  } else if k >= self.n {
    1.0
  } else {
    @special.betainc(
      (self.n - k).to_double(),
      (k + 1).to_double(),
      1.0 - self.p,
    )
  }
}

///|
/// P(X > k) = I_p(k + 1, n - k).
pub impl IntCdf for Binomial with fn sf_int(self, k) {
  if k < 0 {
    1.0
  } else if k >= self.n {
    0.0
  } else {
    @special.betainc((k + 1).to_double(), (self.n - k).to_double(), self.p)
  }
}

///|
/// A normal-approximation start, corrected by exact comparisons.
pub impl IntQuantile for Binomial with fn quantile_int(self, p) {
  let n = self.n.to_double()
  let guess = n * self.p +
    @special.ndtri(p) * (n * self.p * (1.0 - self.p)).sqrt()
  refine_quantile(self, p, @math.floor(guess).to_int(), 0, self.n)
}

///|
pub impl Moments for Binomial with fn mean(self) {
  Some(self.n.to_double() * self.p)
}

///|
pub impl Moments for Binomial with fn variance(self) {
  Some(self.n.to_double() * self.p * (1.0 - self.p))
}

///|
/// φ(t) = (1 - p + p e^{it})^n, raised to the n-th power in polar form.
pub impl Characteristic for Binomial with fn cf(self, t) {
  let re = 1.0 - self.p + self.p * @math.cos(t)
  let im = self.p * @math.sin(t)
  let n = self.n.to_double()
  let r = @math.pow(@math.hypot(re, im), n)
  let theta = n * @math.atan2(im, re)
  (r * @math.cos(theta), r * @math.sin(theta))
}

// ---------------------------------------------------------------------------
// Reference distributions for tests and power
// ---------------------------------------------------------------------------

///|
/// P(X <= x) = P(df/2, x/2), the lower regularised incomplete gamma.
pub impl Cdf for ChiSquared with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @special.gammainc(0.5 * self.df, 0.5 * x)
  }
}

///|
/// P(X > x) = Q(df/2, x/2), computed directly for small p-values.
pub impl Cdf for ChiSquared with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @special.gammaincc(0.5 * self.df, 0.5 * x)
  }
}

///|
pub impl Density for ChiSquared with fn pdf(self, x) {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for ChiSquared with fn logpdf(self, x) {
  if x < 0.0 {
    -1.0 / 0.0
  } else {
    let k = 0.5 * self.df
    (k - 1.0) * @math.ln(x) - 0.5 * x - k * @math.ln(2.0) - @special.lgamma(k)
  }
}

///|
pub impl Quantile for ChiSquared with fn quantile(self, p) {
  2.0 * @special.gammaincinv(0.5 * self.df, p)
}

///|
pub impl Quantile for ChiSquared with fn isf(self, q) {
  2.0 * @special.gammainccinv(0.5 * self.df, q)
}

///|
/// (P(X <= x), P(X > x)) = I_u(dfn/2, dfd/2) and its complement, with
/// u = dfn x/(dfn x + dfd) and 1 - u = dfd/(dfn x + dfd) both formed directly.
fn FisherF::tails(self : FisherF, x : Double) -> (Double, Double) {
  let d = self.dfn * x + self.dfd
  @special.betainc_pair(
    0.5 * self.dfn,
    0.5 * self.dfd,
    self.dfn * x / d,
    self.dfd / d,
  )
}

///|
pub impl Cdf for FisherF with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    self.tails(x).0
  }
}

///|
pub impl Cdf for FisherF with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    self.tails(x).1
  }
}

///|
pub impl Density for FisherF with fn pdf(self, x) {
  if x < 0.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for FisherF with fn logpdf(self, x) {
  if x < 0.0 {
    -1.0 / 0.0
  } else {
    let m = self.dfn
    let n = self.dfd
    0.5 * m * @math.ln(m / n) +
    (0.5 * m - 1.0) * @math.ln(x) -
    0.5 * (m + n) * @math.ln_1p(m * x / n) -
    @special.lbeta(0.5 * m, 0.5 * n)
  }
}

///|
/// The lower tail inverts I_u(dfn/2, dfd/2); above the median the upper
/// tail is inverted instead, so a small 1 - u is never formed.
pub impl Quantile for FisherF with fn quantile(self, p) {
  if p > 0.5 {
    return Quantile::isf(self, 1.0 - p)
  }
  let u = @special.betaincinv(0.5 * self.dfn, 0.5 * self.dfd, p)
  self.dfd * u / (self.dfn * (1.0 - u))
}

///|
pub impl Quantile for FisherF with fn isf(self, q) {
  if q > 0.5 {
    return Quantile::quantile(self, 1.0 - q)
  }
  let v = @special.betaincinv(0.5 * self.dfd, 0.5 * self.dfn, q)
  self.dfd * (1.0 - v) / (self.dfn * v)
}

///|
pub impl Density for Beta with fn pdf(self, x) {
  if x < 0.0 || x > 1.0 {
    0.0
  } else {
    @math.exp(Density::logpdf(self, x))
  }
}

///|
pub impl Density for Beta with fn logpdf(self, x) {
  if x < 0.0 || x > 1.0 {
    -1.0 / 0.0
  } else {
    (self.a - 1.0) * @math.ln(x) +
    (self.b - 1.0) * @math.ln_1p(-x) -
    @special.lbeta(self.a, self.b)
  }
}

///|
pub impl Cdf for Beta with fn sf(self, x) {
  if x <= 0.0 {
    1.0
  } else if x >= 1.0 {
    0.0
  } else {
    @special.betaincc(self.a, self.b, x)
  }
}

///|
/// 1 - I_x(a, b) = q is I_{1-x}(b, a) = q.
pub impl Quantile for Beta with fn isf(self, q) {
  1.0 - @special.betaincinv(self.b, self.a, q)
}

///|
pub impl Cdf for Beta with fn cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else if x >= 1.0 {
    1.0
  } else {
    @special.betainc(self.a, self.b, x)
  }
}

///|
/// Clopper–Pearson intervals for a proportion are quantiles of this.
pub impl Quantile for Beta with fn quantile(self, p) {
  @special.betaincinv(self.a, self.b, p)
}

///|
/// P(X <= k) = Q(k + 1, mu).
pub impl IntCdf for Poisson with fn cdf_int(self, k) {
  if k < 0 {
    0.0
  } else {
    @special.gammaincc((k + 1).to_double(), self.mu)
  }
}

///|
pub impl IntCdf for Poisson with fn sf_int(self, k) {
  if k < 0 {
    1.0
  } else {
    @special.gammainc((k + 1).to_double(), self.mu)
  }
}

///|
pub impl Mass for Poisson with fn pmf(self, k) {
  if k < 0 {
    0.0
  } else if self.mu == 0.0 {
    if k == 0 {
      1.0
    } else {
      0.0
    }
  } else {
    @math.exp(
      k.to_double() * @math.ln(self.mu) -
      self.mu -
      @special.lgamma(k.to_double() + 1.0),
    )
  }
}

///|
pub impl IntQuantile for Poisson with fn quantile_int(self, p) {
  let guess = self.mu + @special.ndtri(p) * self.mu.sqrt()
  refine_quantile(self, p, @math.floor(guess).to_int(), 0, 2147483647)
}

///|
fn Hypergeometric::lowest(self : Hypergeometric) -> Int32 {
  let k = self.draws - (self.population - self.successes)
  if k > 0 {
    k
  } else {
    0
  }
}

///|
fn Hypergeometric::highest(self : Hypergeometric) -> Int32 {
  if self.successes < self.draws {
    self.successes
  } else {
    self.draws
  }
}

///|
/// The null distribution of Fisher's exact test for a 2 × 2 table.
pub impl Mass for Hypergeometric with fn pmf(self, k) {
  if k < self.lowest() || k > self.highest() {
    0.0
  } else {
    @math.exp(
      @special.lchoose(self.successes.to_double(), k.to_double()) +
      @special.lchoose(
        (self.population - self.successes).to_double(),
        (self.draws - k).to_double(),
      ) -
      @special.lchoose(self.population.to_double(), self.draws.to_double()),
    )
  }
}

///|
/// Each tail is summed on its own side, so neither is 1 minus the other.
pub impl IntCdf for Hypergeometric with fn cdf_int(self, k) {
  if k < self.lowest() {
    0.0
  } else if k >= self.highest() {
    1.0
  } else {
    let mut s = @num.CompensatedSum::zero()
    for j in self.lowest()..<=k {
      s = s.add(Mass::pmf(self, j))
    }
    s.value()
  }
}

///|
pub impl IntCdf for Hypergeometric with fn sf_int(self, k) {
  if k < self.lowest() {
    1.0
  } else if k >= self.highest() {
    0.0
  } else {
    let mut s = @num.CompensatedSum::zero()
    for j in (k + 1)..<=self.highest() {
      s = s.add(Mass::pmf(self, j))
    }
    s.value()
  }
}

///|
pub impl IntQuantile for Hypergeometric with fn quantile_int(self, p) {
  let guess = self.draws.to_double() *
    self.successes.to_double() /
    self.population.to_double()
  refine_quantile(
    self,
    p,
    @math.floor(guess).to_int(),
    self.lowest(),
    self.highest(),
  )
}

///|
/// The distribution of a t statistic under the alternative, the basis of
/// t-test power.
pub impl Cdf for NoncentralT with fn cdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl Cdf for NoncentralChiSquared with fn cdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl Cdf for NoncentralF with fn cdf(self, x) {
  ignore((self, x))
  ...
}
