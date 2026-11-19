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
  TooFewPoints(needed~ : Int32, got~ : Int32)
  OutsideSupport(value~ : Double)
  NotConverged(iterations~ : Int32)
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
  fn cdf(Self, Double) -> Double
  fn sf(Self, Double) -> Double = _
  fn logcdf(Self, Double) -> Double = _
  fn logsf(Self, Double) -> Double = _
  /// P(X >= x). Equal to `sf` unless x is an atom; a p-value from a
  /// discrete or randomisation null must count the observed value itself,
  /// so distributions with atoms override this.
  fn sf_at_least(Self, Double) -> Double = _
}

///|
impl Cdf with fn sf_at_least(self, x) {
  self.sf(x)
}

///|
impl Cdf with fn sf(self, x) {
  1.0 - self.cdf(x)
}

///|
impl Cdf with fn logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl Cdf with fn logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
/// The generalised inverse: `quantile(p)` is the smallest x with
/// F(x) >= p, which is defined for atoms and flat stretches too. `isf(q)`
/// inverts `sf` so that an upper-tail probability such as 1e-20 is never
/// rounded through 1 - q.
pub(open) trait Quantile {
  fn quantile(Self, Double) -> Double
  fn isf(Self, Double) -> Double = _
}

///|
impl Quantile with fn isf(self, q) {
  self.quantile(1.0 - q)
}

///|
pub(open) trait Sample {
  fn sample(Self, &@num.RandomSource) -> Double
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
  fn pdf(Self, Double) -> Double
  /// Defaults to log(pdf); families with a closed form override it so that
  /// far tails do not underflow to -inf.
  fn logpdf(Self, Double) -> Double = _
}

///|
impl Density with fn logpdf(self, x) {
  @math.ln(self.pdf(x))
}

///|
/// Point masses on the integers.
pub(open) trait Mass {
  fn pmf(Self, Int32) -> Double
  fn logpmf(Self, Int32) -> Double = _
}

///|
impl Mass with fn logpmf(self, k) {
  @math.ln(self.pmf(k))
}

///|
/// Integer versions of the distribution function and its inverse for
/// distributions on the integers, so callers need not go through Double.
pub(open) trait IntCdf {
  fn cdf_int(Self, Int32) -> Double
  fn sf_int(Self, Int32) -> Double = _
}

///|
impl IntCdf with fn sf_int(self, k) {
  1.0 - self.cdf_int(k)
}

///|
/// The smallest k with P(X <= k) >= p.
pub(open) trait IntQuantile {
  fn quantile_int(Self, Double) -> Int32
}

///|
pub(open) trait SampleInt {
  fn sample_int(Self, &@num.RandomSource) -> Int32
}

// ---------------------------------------------------------------------------
// Other descriptions of a distribution.
// ---------------------------------------------------------------------------

///|
/// `None` when the moment does not exist (Cauchy mean, Student t variance
/// for df <= 1); `Some(inf)` when it exists but is infinite.
pub(open) trait Moments {
  fn mean(Self) -> Double?
  fn variance(Self) -> Double?
}

///|
/// Higher shape descriptors, separate because many families have no
/// closed form for them.
pub(open) trait Shape {
  fn skewness(Self) -> Double?
  fn excess_kurtosis(Self) -> Double?
}

///|
pub(open) trait Entropy {
  /// Differential entropy for densities, Shannon entropy for masses, in nats.
  fn entropy(Self) -> Double
}

///|
/// φ(t) = E[exp(itX)], which exists for every distribution. Returned as
/// (real, imaginary) so this package does not commit to a complex type.
/// Sums of independent variables multiply their characteristic functions,
/// and a distribution can be recovered from φ numerically by an FFT.
pub(open) trait Characteristic {
  fn cf(Self, Double) -> (Double, Double)
}

///|
/// M(t) = E[exp(tX)], `None` where it diverges (it need not exist near any
/// t other than 0, as for the log-normal). The cumulant generating
/// function defaults to its log.
pub(open) trait Mgf {
  fn mgf(Self, Double) -> Double?
  fn cgf(Self, Double) -> Double? = _
}

///|
impl Mgf with fn cgf(self, t) {
  match self.mgf(t) {
    Some(m) => Some(@math.ln(m))
    None => None
  }
}

///|
/// Maximum likelihood (or the family's standard estimator) from a sample.
/// A static method: `Normal::fit(data)` returns a new distribution.
pub(open) trait Fit {
  fn fit(ArrayView[Double]) -> Self raise FitError
}
