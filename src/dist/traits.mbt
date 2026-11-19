///|
/// Capability traits. A distribution implements only what it can compute;
/// an algorithm asks only for what it needs. MoonBit traits have no
/// associated types, so a continuous distribution (support in Double) and a
/// discrete one (support in Int) use separate traits instead of one trait
/// generic over the support type.

///|
/// Raised when a constructor is given parameters outside the family's
/// parameter space. Distributions are immutable, so this is the only place
/// parameters are checked.
#derive(Debug)
pub(all) suberror ParamError {
  InvalidParam(name~ : String, value~ : Double, reason~ : String)
}

///|
/// A source of uniform bits. Sampling takes it as a trait object so that
/// the distribution traits need no type parameter of their own, and so that
/// no global random state exists anywhere in the library.
pub(open) trait RandomSource {
  next_u64(Self) -> UInt64
}

// ---------------------------------------------------------------------------
// Continuous distributions: support in Double.
// ---------------------------------------------------------------------------

///|
pub(open) trait Density {
  pdf(Self, Double) -> Double
  /// Defaults to `pdf(x).ln()`; families with a closed form override it so
  /// that far tails do not underflow to -inf.
  logpdf(Self, Double) -> Double = _
}

///|
impl Density with logpdf(self, x) {
  @math.ln(self.pdf(x))
}

///|
/// Both tails are first class. `sf` defaults to `1 - cdf`, which loses all
/// digits in the upper tail, so every family with a better formula overrides
/// it. The log variants exist for p-values far below 1e-308.
pub(open) trait Cdf {
  cdf(Self, Double) -> Double
  sf(Self, Double) -> Double = _
  logcdf(Self, Double) -> Double = _
  logsf(Self, Double) -> Double = _
}

///|
impl Cdf with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
impl Cdf with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl Cdf with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
/// `quantile(p)` inverts `cdf`, `isf(q)` inverts `sf`. Separate so that an
/// upper-tail probability such as 1e-20 is never rounded through 1 - q.
pub(open) trait Quantile {
  quantile(Self, Double) -> Double
  isf(Self, Double) -> Double = _
}

///|
impl Quantile with isf(self, q) {
  self.quantile(1.0 - q)
}

///|
pub(open) trait SampleContinuous {
  sample(Self, &RandomSource) -> Double
}

// ---------------------------------------------------------------------------
// Discrete distributions: support in Int.
// ---------------------------------------------------------------------------

///|
pub(open) trait Mass {
  pmf(Self, Int) -> Double
  logpmf(Self, Int) -> Double = _
}

///|
impl Mass with logpmf(self, k) {
  @math.ln(self.pmf(k))
}

///|
/// P(X <= k) and P(X > k).
pub(open) trait DiscreteCdf {
  cdf(Self, Int) -> Double
  sf(Self, Int) -> Double = _
}

///|
impl DiscreteCdf with sf(self, k) {
  1.0 - self.cdf(k)
}

///|
/// The smallest k with cdf(k) >= p.
pub(open) trait DiscreteQuantile {
  quantile(Self, Double) -> Int
}

///|
pub(open) trait SampleDiscrete {
  sample(Self, &RandomSource) -> Int
}

// ---------------------------------------------------------------------------
// Shared.
// ---------------------------------------------------------------------------

///|
/// `None` when the moment does not exist (Cauchy mean, Student t variance
/// for df <= 1); `Some(inf)` when it exists but is infinite.
pub(open) trait Moments {
  mean(Self) -> Double?
  variance(Self) -> Double?
}
