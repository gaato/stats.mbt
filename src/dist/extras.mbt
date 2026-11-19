///|
/// Capabilities written by hand on top of the generated catalogue types:
/// the special-function paths of the reference distributions, moment and
/// characteristic functions, fitting, and alternative constructors. Field
/// names follow SciPy (generated in catalog.mbt).

// ---------------------------------------------------------------------------
// Normal: loc and scale built in.
// ---------------------------------------------------------------------------

///|
pub impl Density for Normal with pdf(self, x) {
  @math.exp(Density::logpdf(self, x))
}

///|
/// Closed form, so the far tails stay finite on the log scale.
pub impl Density for Normal with logpdf(self, x) {
  let z = (x - self.loc) / self.scale
  -0.5 * z * z - @math.ln(self.scale) - 0.5 * @math.ln(2.0 * @math.PI)
}

///|
pub impl Cdf for Normal with cdf(self, x) {
  @special.ndtr((x - self.loc) / self.scale)
}

///|
/// Symmetry gives the upper tail without cancellation.
pub impl Cdf for Normal with sf(self, x) {
  @special.ndtr((self.loc - x) / self.scale)
}

///|
pub impl Quantile for Normal with quantile(self, p) {
  self.loc + self.scale * @special.ndtri(p)
}

///|
pub impl Moments for Normal with mean(self) {
  Some(self.loc)
}

///|
pub impl Moments for Normal with variance(self) {
  Some(self.scale * self.scale)
}

///|
/// φ(t) = exp(i loc t - scale² t²/2).
pub impl Characteristic for Normal with cf(self, t) {
  let r = @math.exp(-0.5 * self.scale * self.scale * t * t)
  (r * @math.cos(self.loc * t), r * @math.sin(self.loc * t))
}

///|
/// M(t) = exp(loc t + scale² t²/2), finite for every t.
pub impl Mgf for Normal with mgf(self, t) {
  Some(@math.exp(self.loc * t + 0.5 * self.scale * self.scale * t * t))
}

///|
/// Maximum likelihood: the sample mean and the divisor-n standard
/// deviation. Needs two points, since one gives scale = 0.
pub impl Fit for Normal with fit(data) {
  if data.length() < 2 {
    raise TooFewPoints(needed=2, got=data.length())
  }
  ...
}

// ---------------------------------------------------------------------------
// Student t
// ---------------------------------------------------------------------------

///|
pub impl Moments for StudentT with mean(self) {
  if self.df > 1.0 {
    Some(0.0)
  } else {
    None
  }
}

///|
pub impl Moments for StudentT with variance(self) {
  if self.df > 2.0 {
    Some(self.df / (self.df - 2.0))
  } else if self.df > 1.0 {
    Some(1.0 / 0.0)
  } else {
    None
  }
}

///|
/// No `Mgf`: E[exp(tX)] diverges for every t != 0. The characteristic
/// function is a modified Bessel function of the second kind.
pub impl Characteristic for StudentT with cf(self, t) {
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
      { name: "a", constraint: Positive },
      { name: "rate", constraint: Positive },
    ],
    make: p => Gamma::with_rate(p[0], p[1]),
    params: d => [d.a, 1.0 / d.scale],
  }
}

///|
pub impl Cdf for Gamma with cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @special.gammainc(self.a, x / self.scale)
  }
}

///|
pub impl Cdf for Gamma with sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @special.gammaincc(self.a, x / self.scale)
  }
}

///|
/// M(t) = (1 - t scale)^(-a) for t < 1/scale; it diverges from there on.
pub impl Mgf for Gamma with mgf(self, t) {
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
pub impl Mass for Binomial with pmf(self, k) {
  ignore((self, k))
  ...
}

///|
/// P(X <= k) = I_{1-p}(n - k, k + 1).
pub impl IntCdf for Binomial with cdf_int(self, k) {
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
pub impl IntCdf for Binomial with sf_int(self, k) {
  if k < 0 {
    1.0
  } else if k >= self.n {
    0.0
  } else {
    @special.betainc((k + 1).to_double(), (self.n - k).to_double(), self.p)
  }
}

///|
pub impl IntQuantile for Binomial with quantile_int(self, p) {
  ignore((self, p))
  ...
}

///|
pub impl Moments for Binomial with mean(self) {
  Some(self.n.to_double() * self.p)
}

///|
pub impl Moments for Binomial with variance(self) {
  Some(self.n.to_double() * self.p * (1.0 - self.p))
}

///|
/// φ(t) = (1 - p + p e^{it})^n, raised to the n-th power in polar form.
pub impl Characteristic for Binomial with cf(self, t) {
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
pub impl Cdf for ChiSquared with cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @special.gammainc(0.5 * self.df, 0.5 * x)
  }
}

///|
/// P(X > x) = Q(df/2, x/2), computed directly for small p-values.
pub impl Cdf for ChiSquared with sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @special.gammaincc(0.5 * self.df, 0.5 * x)
  }
}

///|
/// P(X > x) = I_{dfd/(dfd + dfn x)}(dfd/2, dfn/2).
pub impl Cdf for FisherF with sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @special.betainc(
      0.5 * self.dfd,
      0.5 * self.dfn,
      self.dfd / (self.dfd + self.dfn * x),
    )
  }
}

///|
pub impl Cdf for FisherF with cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @special.betaincc(
      0.5 * self.dfd,
      0.5 * self.dfn,
      self.dfd / (self.dfd + self.dfn * x),
    )
  }
}

///|
pub impl Cdf for Beta with cdf(self, x) {
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
pub impl Quantile for Beta with quantile(self, p) {
  @special.betaincinv(self.a, self.b, p)
}

///|
/// P(X <= k) = Q(k + 1, mu).
pub impl IntCdf for Poisson with cdf_int(self, k) {
  if k < 0 {
    0.0
  } else {
    @special.gammaincc((k + 1).to_double(), self.mu)
  }
}

///|
pub impl IntCdf for Poisson with sf_int(self, k) {
  if k < 0 {
    1.0
  } else {
    @special.gammainc((k + 1).to_double(), self.mu)
  }
}

///|
/// The null distribution of Fisher's exact test for a 2 × 2 table.
pub impl Mass for Hypergeometric with pmf(self, k) {
  ignore((self, k))
  ...
}

///|
/// The distribution of a t statistic under the alternative, the basis of
/// t-test power.
pub impl Cdf for NoncentralT with cdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl Cdf for NoncentralChiSquared with cdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl Cdf for NoncentralF with cdf(self, x) {
  ignore((self, x))
  ...
}
