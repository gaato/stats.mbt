///|
/// Families. Each is an immutable struct with private fields; constructors
/// check the parameter space and raise `ParamError`. Where conventions
/// differ between libraries (scale vs rate), the constructor name says
/// which one is meant, and there is no constructor that guesses.

// ---------------------------------------------------------------------------
// Normal
// ---------------------------------------------------------------------------

///|
pub struct Normal {
  priv mu : Double
  priv sigma : Double
}

///|
pub fn Normal::new(
  mu? : Double = 0.0,
  sigma? : Double = 1.0,
) -> Normal raise ParamError {
  if !(sigma > 0.0) {
    raise InvalidParam(name="sigma", value=sigma, reason="must be > 0")
  }
  { mu, sigma }
}

///|
pub impl Density for Normal with pdf(self, x) {
  @math.exp(Density::logpdf(self, x))
}

///|
/// Closed form, so the far tails stay finite on the log scale.
pub impl Density for Normal with logpdf(self, x) {
  let z = (x - self.mu) / self.sigma
  -0.5 * z * z - @math.ln(self.sigma) - 0.5 * @math.ln(2.0 * @math.PI)
}

///|
pub impl Cdf for Normal with cdf(self, x) {
  @special.ndtr((x - self.mu) / self.sigma)
}

///|
/// Overridden: symmetry gives the upper tail without cancellation.
pub impl Cdf for Normal with sf(self, x) {
  @special.ndtr((self.mu - x) / self.sigma)
}

///|
pub impl Quantile for Normal with quantile(self, p) {
  self.mu + self.sigma * @special.ndtri(p)
}

///|
pub impl Sample for Normal with sample(self, rng) {
  ignore((self, rng))
  ...
}

///|
pub impl Univariate for Normal

///|
pub impl Moments for Normal with mean(self) {
  Some(self.mu)
}

///|
pub impl Moments for Normal with variance(self) {
  Some(self.sigma * self.sigma)
}

///|
/// φ(t) = exp(iμt - σ²t²/2).
pub impl Characteristic for Normal with cf(self, t) {
  let r = @math.exp(-0.5 * self.sigma * self.sigma * t * t)
  (r * @math.cos(self.mu * t), r * @math.sin(self.mu * t))
}

///|
/// M(t) = exp(μt + σ²t²/2), finite for every t.
pub impl Mgf for Normal with mgf(self, t) {
  Some(@math.exp(self.mu * t + 0.5 * self.sigma * self.sigma * t * t))
}

///|
/// Maximum likelihood: the sample mean and the divisor-n standard
/// deviation. Needs two points, since one gives sigma = 0.
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
pub struct StudentT {
  priv df : Double
}

///|
pub fn StudentT::new(df : Double) -> StudentT raise ParamError {
  if !(df > 0.0) {
    raise InvalidParam(name="df", value=df, reason="must be > 0")
  }
  { df, }
}

///|
pub impl Density for StudentT with pdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl Cdf for StudentT with cdf(self, x) {
  ignore((self, x))
  ...
}

///|
/// The p-value path of every t test. Computed from the incomplete beta
/// directly, never as 1 - cdf.
pub impl Cdf for StudentT with sf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl Quantile for StudentT with quantile(self, p) {
  ignore((self, p))
  ...
}

///|
pub impl Sample for StudentT with sample(self, rng) {
  ignore((self, rng))
  ...
}

///|
pub impl Univariate for StudentT

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
/// function exists (a modified Bessel function of the second kind), so
/// `Characteristic` can come later without changing anything else.
pub impl Characteristic for StudentT with cf(self, t) {
  ignore((self, t))
  ...
}

// ---------------------------------------------------------------------------
// Gamma
// ---------------------------------------------------------------------------

///|
/// Stored as shape and rate; both constructors exist and neither is the
/// default, because SciPy (scale) and much of the literature (rate) disagree.
pub struct Gamma {
  priv shape : Double
  priv rate : Double
}

///|
pub fn Gamma::with_rate(
  shape : Double,
  rate : Double,
) -> Gamma raise ParamError {
  if !(shape > 0.0) {
    raise InvalidParam(name="shape", value=shape, reason="must be > 0")
  }
  if !(rate > 0.0) {
    raise InvalidParam(name="rate", value=rate, reason="must be > 0")
  }
  { shape, rate }
}

///|
pub fn Gamma::with_scale(
  shape : Double,
  scale : Double,
) -> Gamma raise ParamError {
  if !(scale > 0.0) {
    raise InvalidParam(name="scale", value=scale, reason="must be > 0")
  }
  Gamma::with_rate(shape, 1.0 / scale)
}

///|
pub impl Cdf for Gamma with cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @special.gammainc(self.shape, self.rate * x)
  }
}

///|
pub impl Cdf for Gamma with sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @special.gammaincc(self.shape, self.rate * x)
  }
}

///|
/// M(t) = (1 - t/rate)^(-shape) for t < rate; it diverges from t = rate on.
pub impl Mgf for Gamma with mgf(self, t) {
  if t < self.rate {
    Some(@math.pow(1.0 - t / self.rate, -self.shape))
  } else {
    None
  }
}

// ---------------------------------------------------------------------------
// Binomial
// ---------------------------------------------------------------------------

///|
pub struct Binomial {
  priv n : Int
  priv p : Double
}

///|
pub fn Binomial::new(n : Int, p : Double) -> Binomial raise ParamError {
  if n < 0 {
    raise InvalidParam(name="n", value=n.to_double(), reason="must be >= 0")
  }
  if !(p >= 0.0 && p <= 1.0) {
    raise InvalidParam(name="p", value=p, reason="must be in [0, 1]")
  }
  { n, p }
}

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
pub impl SampleInt for Binomial with sample_int(self, rng) {
  ignore((self, rng))
  ...
}

///|
/// On the real line the distribution function is the step function.
pub impl Cdf for Binomial with cdf(self, x) {
  IntCdf::cdf_int(self, @math.floor(x).to_int())
}

///|
pub impl Cdf for Binomial with sf(self, x) {
  IntCdf::sf_int(self, @math.floor(x).to_int())
}

///|
/// P(X >= x) = P(X > ceil(x) - 1).
pub impl Cdf for Binomial with sf_at_least(self, x) {
  IntCdf::sf_int(self, @math.ceil(x).to_int() - 1)
}

///|
pub impl Quantile for Binomial with quantile(self, p) {
  IntQuantile::quantile_int(self, p).to_double()
}

///|
pub impl Sample for Binomial with sample(self, rng) {
  SampleInt::sample_int(self, rng).to_double()
}

///|
pub impl Univariate for Binomial

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
