///|
/// Capability traits for distributions on the real line.
///
/// Every real-valued distribution has a distribution function, whatever its
/// Lebesgue decomposition (atoms, a density, a singular part). So `Cdf`,
/// `Quantile` and `Sample` on `Double` are the common capabilities, and a
/// density or a mass function is an extra capability that only some
/// distributions have. A distribution on the integers also implements the
/// `Int` traits for convenience; its `Double` distribution function is the
/// step function.

///|
/// Raised when a constructor is given parameters outside the family's
/// parameter space. Distributions are immutable, so this is the only place
/// parameters are checked.
#derive(Debug)
pub(all) suberror ParamError {
  InvalidParam(name~ : String, value~ : Double, reason~ : String)
}

///|
/// Raised when fitting a family to data fails: too few points, data outside
/// the support, or an optimiser that did not converge. A fit never returns
/// parameters it could not justify.
#derive(Debug)
pub(all) suberror FitError {
  TooFewPoints(needed~ : Int, got~ : Int)
  OutsideSupport(value~ : Double)
  NotConverged(iterations~ : Int)
}

///|
/// A source of uniform bits. Sampling takes it as a trait object so that the
/// distribution traits need no type parameter of their own, and so that no
/// global random state exists anywhere in the library.
pub(open) trait RandomSource {
  next_u64(Self) -> UInt64
}

// ---------------------------------------------------------------------------
// Common to every distribution on the real line.
// ---------------------------------------------------------------------------

///|
/// F(x) = P(X <= x). Both tails are first class: `sf` defaults to `1 - cdf`,
/// which loses every digit in the upper tail, so a family with a better
/// formula overrides it. The log variants exist for p-values below the
/// smallest positive Double.
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
/// The generalised inverse: `quantile(p)` is the smallest x with
/// F(x) >= p, which is defined for atoms and flat stretches too. `isf(q)`
/// inverts `sf` so that an upper-tail probability such as 1e-20 is never
/// rounded through 1 - q.
pub(open) trait Quantile {
  quantile(Self, Double) -> Double
  isf(Self, Double) -> Double = _
}

///|
impl Quantile with isf(self, q) {
  self.quantile(1.0 - q)
}

///|
pub(open) trait Sample {
  sample(Self, &RandomSource) -> Double
}

///|
/// The bundle a composite (a mixture, a truncation) needs from its parts.
/// Used as a trait object `&Univariate`; through the object the methods are
/// called as `Cdf::cdf(d, x)`, since dot syntax does not reach supertraits.
pub(open) trait Univariate: Cdf + Quantile + Sample {}

// ---------------------------------------------------------------------------
// Parts of the Lebesgue decomposition that only some distributions have.
// ---------------------------------------------------------------------------

///|
/// The absolutely continuous part has a density.
pub(open) trait Density {
  pdf(Self, Double) -> Double
  /// Defaults to log(pdf); families with a closed form override it so that
  /// far tails do not underflow to -inf.
  logpdf(Self, Double) -> Double = _
}

///|
impl Density with logpdf(self, x) {
  @math.ln(self.pdf(x))
}

///|
/// Point masses on the integers.
pub(open) trait Mass {
  pmf(Self, Int) -> Double
  logpmf(Self, Int) -> Double = _
}

///|
impl Mass with logpmf(self, k) {
  @math.ln(self.pmf(k))
}

///|
/// Integer versions of the distribution function and its inverse for
/// distributions on the integers, so callers need not go through Double.
pub(open) trait IntCdf {
  cdf_int(Self, Int) -> Double
  sf_int(Self, Int) -> Double = _
}

///|
impl IntCdf with sf_int(self, k) {
  1.0 - self.cdf_int(k)
}

///|
/// The smallest k with P(X <= k) >= p.
pub(open) trait IntQuantile {
  quantile_int(Self, Double) -> Int
}

///|
pub(open) trait SampleInt {
  sample_int(Self, &RandomSource) -> Int
}

// ---------------------------------------------------------------------------
// Other descriptions of a distribution.
// ---------------------------------------------------------------------------

///|
/// `None` when the moment does not exist (Cauchy mean, Student t variance
/// for df <= 1); `Some(inf)` when it exists but is infinite.
pub(open) trait Moments {
  mean(Self) -> Double?
  variance(Self) -> Double?
}

///|
/// Higher shape descriptors, separate because many families have no
/// closed form for them.
pub(open) trait Shape {
  skewness(Self) -> Double?
  excess_kurtosis(Self) -> Double?
}

///|
pub(open) trait Entropy {
  /// Differential entropy for densities, Shannon entropy for masses, in nats.
  entropy(Self) -> Double
}

///|
/// φ(t) = E[exp(itX)], which exists for every distribution. Returned as
/// (real, imaginary) so this package does not commit to a complex type.
/// Sums of independent variables multiply their characteristic functions,
/// and a distribution can be recovered from φ numerically by an FFT.
pub(open) trait Characteristic {
  cf(Self, Double) -> (Double, Double)
}

///|
/// M(t) = E[exp(tX)], `None` where it diverges (it need not exist near any
/// t other than 0, as for the log-normal). The cumulant generating
/// function defaults to its log.
pub(open) trait Mgf {
  mgf(Self, Double) -> Double?
  cgf(Self, Double) -> Double? = _
}

///|
impl Mgf with cgf(self, t) {
  match self.mgf(t) {
    Some(m) => Some(@math.ln(m))
    None => None
  }
}

///|
/// Maximum likelihood (or the family's standard estimator) from a sample.
/// A static method: `Normal::fit(data)` returns a new distribution.
pub(open) trait Fit {
  fit(ArrayView[Double]) -> Self raise FitError
}
