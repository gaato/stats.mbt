///|
/// Families. Each is an immutable struct with private fields; constructors
/// check the parameter space and raise `ParamError`. Where conventions
/// differ between libraries (scale vs rate), the constructor name says
/// which one is meant, and there is no constructor that guesses.

///|
pub struct Normal {
  priv mu : Double
  priv sigma : Double
}

///|
pub fn Normal::new(mu? : Double = 0.0, sigma? : Double = 1.0) -> Normal raise ParamError {
  if !(sigma > 0.0) {
    raise InvalidParam(name="sigma", value=sigma, reason="must be > 0")
  }
  { mu, sigma }
}

///|
pub impl Density for Normal with pdf(self, x) {
  ignore((self, x))
  ...
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
pub impl Moments for Normal with mean(self) {
  Some(self.mu)
}

///|
pub impl Moments for Normal with variance(self) {
  Some(self.sigma * self.sigma)
}

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
/// Stored as shape and rate; both constructors exist and neither is the
/// default, because SciPy (scale) and much of the literature (rate) disagree.
pub struct Gamma {
  priv shape : Double
  priv rate : Double
}

///|
pub fn Gamma::with_rate(shape : Double, rate : Double) -> Gamma raise ParamError {
  if !(shape > 0.0) {
    raise InvalidParam(name="shape", value=shape, reason="must be > 0")
  }
  if !(rate > 0.0) {
    raise InvalidParam(name="rate", value=rate, reason="must be > 0")
  }
  { shape, rate }
}

///|
pub fn Gamma::with_scale(shape : Double, scale : Double) -> Gamma raise ParamError {
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
pub impl DiscreteCdf for Binomial with cdf(self, k) {
  if k < 0 {
    0.0
  } else if k >= self.n {
    1.0
  } else {
    // P(X <= k) = I_{1-p}(n - k, k + 1)
    @special.betainc((self.n - k).to_double(), (k + 1).to_double(), 1.0 - self.p)
  }
}

///|
pub impl DiscreteCdf for Binomial with sf(self, k) {
  if k < 0 {
    1.0
  } else if k >= self.n {
    0.0
  } else {
    @special.betainc((k + 1).to_double(), (self.n - k).to_double(), self.p)
  }
}
