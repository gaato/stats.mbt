///|
/// The reference distributions the tests and power analysis need. Bodies
/// come with the special functions; the constructors check parameters now.

///|
pub struct ChiSquared {
  priv df : Double
}

///|
pub fn ChiSquared::new(df : Double) -> ChiSquared raise ParamError {
  if !(df > 0.0) {
    raise InvalidParam(name="df", value=df, reason="must be > 0")
  }
  { df, }
}

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
pub struct FisherF {
  priv df1 : Double
  priv df2 : Double
}

///|
pub fn FisherF::new(df1 : Double, df2 : Double) -> FisherF raise ParamError {
  if !(df1 > 0.0) {
    raise InvalidParam(name="df1", value=df1, reason="must be > 0")
  }
  if !(df2 > 0.0) {
    raise InvalidParam(name="df2", value=df2, reason="must be > 0")
  }
  { df1, df2 }
}

///|
/// P(X > x) = I_{df2/(df2 + df1 x)}(df2/2, df1/2).
pub impl Cdf for FisherF with sf(self, x) {
  if x <= 0.0 {
    1.0
  } else {
    @special.betainc(
      0.5 * self.df2,
      0.5 * self.df1,
      self.df2 / (self.df2 + self.df1 * x),
    )
  }
}

///|
pub impl Cdf for FisherF with cdf(self, x) {
  if x <= 0.0 {
    0.0
  } else {
    @special.betaincc(
      0.5 * self.df2,
      0.5 * self.df1,
      self.df2 / (self.df2 + self.df1 * x),
    )
  }
}

///|
pub struct Beta {
  priv a : Double
  priv b : Double
}

///|
pub fn Beta::new(a : Double, b : Double) -> Beta raise ParamError {
  if !(a > 0.0) {
    raise InvalidParam(name="a", value=a, reason="must be > 0")
  }
  if !(b > 0.0) {
    raise InvalidParam(name="b", value=b, reason="must be > 0")
  }
  { a, b }
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
pub struct Poisson {
  priv lambda : Double
}

///|
pub fn Poisson::new(lambda : Double) -> Poisson raise ParamError {
  if !(lambda >= 0.0) || lambda.is_inf() {
    raise InvalidParam(name="lambda", value=lambda, reason="must be >= 0")
  }
  { lambda, }
}

///|
/// P(X <= k) = Q(k + 1, λ).
pub impl IntCdf for Poisson with cdf_int(self, k) {
  if k < 0 {
    0.0
  } else {
    @special.gammaincc((k + 1).to_double(), self.lambda)
  }
}

///|
pub impl IntCdf for Poisson with sf_int(self, k) {
  if k < 0 {
    1.0
  } else {
    @special.gammainc((k + 1).to_double(), self.lambda)
  }
}

///|
/// Draws of n from a population of N containing K successes: the null of
/// Fisher's exact test for a 2 × 2 table.
pub struct Hypergeometric {
  priv population : Int
  priv successes : Int
  priv draws : Int
}

///|
pub fn Hypergeometric::new(
  population~ : Int,
  successes~ : Int,
  draws~ : Int,
) -> Hypergeometric raise ParamError {
  if population < 0 || successes < 0 || successes > population ||
    draws < 0 ||
    draws > population {
    raise InvalidParam(
      name="successes",
      value=successes.to_double(),
      reason="need 0 <= successes, draws <= population",
    )
  }
  { population, successes, draws }
}

///|
pub impl Mass for Hypergeometric with pmf(self, k) {
  ignore((self, k))
  ...
}

///|
/// Noncentral t: the distribution of a t statistic under the alternative,
/// the basis of t-test power.
pub struct NoncentralT {
  priv df : Double
  priv delta : Double
}

///|
pub fn NoncentralT::new(df : Double, delta : Double) -> NoncentralT raise ParamError {
  if !(df > 0.0) {
    raise InvalidParam(name="df", value=df, reason="must be > 0")
  }
  { df, delta }
}

///|
pub impl Cdf for NoncentralT with cdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub struct NoncentralChiSquared {
  priv df : Double
  priv lambda : Double
}

///|
pub fn NoncentralChiSquared::new(
  df : Double,
  lambda : Double,
) -> NoncentralChiSquared raise ParamError {
  if !(df > 0.0) {
    raise InvalidParam(name="df", value=df, reason="must be > 0")
  }
  if !(lambda >= 0.0) {
    raise InvalidParam(name="lambda", value=lambda, reason="must be >= 0")
  }
  { df, lambda }
}

///|
pub impl Cdf for NoncentralChiSquared with cdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub struct NoncentralF {
  priv df1 : Double
  priv df2 : Double
  priv lambda : Double
}

///|
pub fn NoncentralF::new(
  df1 : Double,
  df2 : Double,
  lambda : Double,
) -> NoncentralF raise ParamError {
  if !(df1 > 0.0) || !(df2 > 0.0) || !(lambda >= 0.0) {
    raise InvalidParam(name="df1", value=df1, reason="need df1, df2 > 0, lambda >= 0")
  }
  { df1, df2, lambda }
}

///|
pub impl Cdf for NoncentralF with cdf(self, x) {
  ignore((self, x))
  ...
}
