///|
/// Distributions on R^d and on Z^d. The univariate distribution function
/// is not the central description here (the joint CDF is rarely what
/// anyone wants), so the common capabilities are the dimension, a density
/// or mass, sampling and the first two moments. Moments come back as
/// values, like the univariate `Moments`; sampling can also write into a
/// caller's buffer, so that a million points need no allocation per point.
/// Matrices are `@data.DataMatrix` until a dense linear algebra module
/// provides a matrix type.

///|
/// Multivariate distributions over vectors of a fixed length.
pub(open) trait Dimension {
  fn dim(Self) -> Int32
}

///|
/// Multivariate distributions with a density.
pub(open) trait MvDensity: Dimension {
  fn logpdf(Self, ArrayView[Double]) -> Double
}

///|
/// Multivariate distributions over vectors of integers.
pub(open) trait MvMass: Dimension {
  fn logpmf(Self, ArrayView[Int64]) -> Double
}

///|
/// One draw. `sample_into` writes it into `out`, which must have length
/// `dim`; `sample` allocates the vector.
pub(open) trait MvSample: Dimension {
  fn sample_into(Self, &@num.RandomSource, FixedArray[Double]) -> Unit
  fn sample(Self, &@num.RandomSource) -> FixedArray[Double] = _
}

///|
impl MvSample with fn sample(self, rng) {
  let out = FixedArray::make(self.dim(), 0.0)
  self.sample_into(rng, out)
  out
}

///|
/// The mean vector (length d) and the covariance matrix (d × d); `None`
/// when the moment does not exist, as for the univariate `Moments`.
pub(open) trait MvMoments: Dimension {
  fn mean(Self) -> FixedArray[Double]?
  fn covariance(Self) -> @data.DataMatrix?
}

///|
/// The univariate marginal of one coordinate, as a trait object so that
/// families with different marginal types share the capability.
pub(open) trait Marginal: Dimension {
  fn marginal(Self, Int32) -> &@dist.Univariate
}

///|
/// Multivariate normal. Stored as the mean and the lower Cholesky factor
/// of the covariance, so that the density and sampling never form an
/// inverse. The constructor factorises, and fails if the covariance is not
/// positive definite.
pub struct MvNormal {
  priv d : Int32
  priv mean : FixedArray[Double]
  priv chol : FixedArray[Double]
}

///|
/// The multivariate normal with this mean and a symmetric positive
/// definite covariance of matching size.
pub fn MvNormal::new(
  mean : ArrayView[Double],
  covariance : @data.DataMatrix,
) -> MvNormal raise @dist.ParamError {
  let d = mean.length()
  if covariance.rows != d || covariance.cols != d {
    raise @dist.InvalidParam(
      name="covariance",
      value=covariance.rows.to_double(),
      reason="must be d x d for a mean of length d",
    )
  }
  ...
}

///|
/// The length of the vectors.
pub impl Dimension for MvNormal with fn dim(self) {
  self.d
}

///|
/// The log density at x.
pub impl MvDensity for MvNormal with fn logpdf(self, x) {
  ignore((self, x))
  ...
}

///|
/// Writes one draw from `rng` into `out`.
pub impl MvSample for MvNormal with fn sample_into(self, rng, out) {
  ignore((self, rng, out))
  ...
}

///|
/// The mean, or `None` when it does not exist.
pub impl MvMoments for MvNormal with fn mean(self) {
  Some(self.mean.copy())
}

///|
/// L Lᵀ from the stored Cholesky factor.
pub impl MvMoments for MvNormal with fn covariance(self) {
  ignore(self)
  ...
}

///|
/// The Dirichlet distribution on the probability simplex.
pub struct Dirichlet {
  priv alpha : FixedArray[Double]
}

///|
/// The Dirichlet distribution with positive concentrations `alpha`.
pub fn Dirichlet::new(
  alpha : ArrayView[Double],
) -> Dirichlet raise @dist.ParamError {
  for a in alpha {
    if !(a > 0.0) {
      raise @dist.InvalidParam(name="alpha", value=a, reason="must be > 0")
    }
  }
  { alpha: FixedArray::from_iter(alpha.iter()), }
}

///|
/// The length of the vectors.
pub impl Dimension for Dirichlet with fn dim(self) {
  self.alpha.length()
}

///|
/// The log density at x.
pub impl MvDensity for Dirichlet with fn logpdf(self, x) {
  ignore((self, x))
  ...
}

///|
/// Writes one draw from `rng` into `out`.
pub impl MvSample for Dirichlet with fn sample_into(self, rng, out) {
  ignore((self, rng, out))
  ...
}

///|
/// Counts of n trials over d categories.
pub struct Multinomial {
  priv n : Int64
  priv p : FixedArray[Double]
}

///|
/// The multinomial distribution of `n` trials over categories with
/// probabilities `p`.
pub fn Multinomial::new(
  n : Int64,
  p : ArrayView[Double],
) -> Multinomial raise @dist.ParamError {
  if n < 0 {
    raise @dist.InvalidParam(
      name="n",
      value=n.to_double(),
      reason="must be >= 0",
    )
  }
  ignore(p)
  ...
}

///|
/// The number of categories.
pub impl Dimension for Multinomial with fn dim(self) {
  self.p.length()
}

///|
/// ln P(X = k) for a vector of counts.
pub impl MvMass for Multinomial with fn logpmf(self, k) {
  ignore((self, k))
  ...
}

///|
/// A copula: a distribution on [0, 1]^d with uniform marginals. Combined
/// with any univariate marginals it builds a joint distribution whose
/// dependence is separate from its marginals (Sklar's theorem).
pub(open) trait Copula: Dimension {
  fn copula_cdf(Self, ArrayView[Double]) -> Double
  fn sample_uniform_into(Self, &@num.RandomSource, FixedArray[Double]) -> Unit
}

///|
/// The Gaussian copula with a given correlation matrix.
pub struct GaussianCopula {
  priv mvn : MvNormal
}

///|
/// The Gaussian copula with a correlation matrix.
pub fn GaussianCopula::new(
  correlation : @data.DataMatrix,
) -> GaussianCopula raise @dist.ParamError {
  ignore(correlation)
  ...
}

///|
/// The number of variables.
pub impl Dimension for GaussianCopula with fn dim(self) {
  self.mvn.d
}

///|
/// C(u), the joint CDF of the uniform margins.
pub impl Copula for GaussianCopula with fn copula_cdf(self, u) {
  ignore((self, u))
  ...
}

///|
/// Writes one draw of uniform margins from `rng` into `out`.
pub impl Copula for GaussianCopula with fn sample_uniform_into(self, rng, out) {
  ignore((self, rng, out))
  ...
}

///|
/// A joint distribution from a copula and one marginal per coordinate.
pub struct Joint[C] {
  priv copula : C
  priv marginals : Array[&@dist.Univariate]
}

///|
/// The joint distribution of a copula and one marginal per variable.
pub fn[C : Copula] Joint::new(
  copula : C,
  marginals : Array[&@dist.Univariate],
) -> Joint[C] raise @dist.ParamError {
  if Dimension::dim(copula) != marginals.length() {
    raise @dist.InvalidParam(
      name="marginals",
      value=marginals.length().to_double(),
      reason="need one marginal per copula dimension",
    )
  }
  { copula, marginals, }
}

///|
/// The number of variables.
pub impl[C : Copula] Dimension for Joint[C] with fn dim(self) {
  self.marginals.length()
}

///|
/// Sample uniforms from the copula, then map each through its marginal's
/// quantile function.
pub impl[C : Copula] MvSample for Joint[C] with fn sample_into(self, rng, out) {
  self.copula.sample_uniform_into(rng, out)
  for i, m in self.marginals {
    out[i] = @dist.Quantile::quantile(m, out[i])
  }
}

///|
/// Dot syntax for `Dimension`: `d.dim`.
pub extend MvNormal with Dimension::{dim}

///|
/// Dot syntax for `MvDensity`: `d.logpdf`.
pub extend MvNormal with MvDensity::{logpdf}

///|
/// Dot syntax for `MvSample`: `d.sample_into`, `d.sample`.
pub extend MvNormal with MvSample::{sample_into, sample}

///|
/// Dot syntax for `MvMoments`: `d.mean`, `d.covariance`.
pub extend MvNormal with MvMoments::{mean, covariance}

///|
/// Dot syntax for `Dimension`: `d.dim`.
pub extend Dirichlet with Dimension::{dim}

///|
/// Dot syntax for `MvDensity`: `d.logpdf`.
pub extend Dirichlet with MvDensity::{logpdf}

///|
/// Dot syntax for `MvSample`: `d.sample_into`, `d.sample`.
pub extend Dirichlet with MvSample::{sample_into, sample}

///|
/// Dot syntax for `Dimension`: `d.dim`.
pub extend Multinomial with Dimension::{dim}

///|
/// Dot syntax for `MvMass`: `d.logpmf`.
pub extend Multinomial with MvMass::{logpmf}

///|
/// Dot syntax for `Dimension`: `d.dim`.
pub extend GaussianCopula with Dimension::{dim}

///|
/// Dot syntax for `Copula`: `d.copula_cdf`, `d.sample_uniform_into`.
pub extend GaussianCopula with Copula::{copula_cdf, sample_uniform_into}

///|
/// Dot syntax for `Dimension`: `d.dim`.
pub extend Joint with Dimension::{dim}

///|
/// Dot syntax for `MvSample`: `d.sample_into`, `d.sample`.
pub extend Joint with MvSample::{sample_into, sample}
