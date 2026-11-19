///|
/// Distributions on R^d and on Z^d. The univariate distribution function
/// is not the central description here (the joint CDF is rarely what
/// anyone wants), so the common capabilities are the dimension, a density
/// or mass, sampling and the first two moments. Vectors are passed as
/// views and written into caller-provided buffers, so sampling a million
/// points allocates nothing per point. Matrices are row-major
/// `FixedArray[Double]` of size d * d until the dense linear algebra
/// module provides a matrix type.

///|
pub(open) trait Dimension {
  dim(Self) -> Int
}

///|
pub(open) trait MvDensity: Dimension {
  logpdf(Self, ArrayView[Double]) -> Double
}

///|
pub(open) trait MvMass: Dimension {
  logpmf(Self, ArrayView[Int]) -> Double
}

///|
/// Writes one draw into `out`, which must have length `dim`.
pub(open) trait MvSample: Dimension {
  sample_into(Self, &@num.RandomSource, FixedArray[Double]) -> Unit
}

///|
/// Writes the mean vector (length d) and the covariance matrix (d * d,
/// row-major) into caller buffers. `false` when a moment does not exist.
pub(open) trait MvMoments: Dimension {
  mean_into(Self, FixedArray[Double]) -> Bool
  covariance_into(Self, FixedArray[Double]) -> Bool
}

///|
/// The univariate marginal of one coordinate, as a trait object so that
/// families with different marginal types share the capability.
pub(open) trait Marginal: Dimension {
  marginal(Self, Int) -> &@dist.Univariate
}

///|
/// Multivariate normal. Stored as the mean and the lower Cholesky factor
/// of the covariance, so that the density and sampling never form an
/// inverse. The constructor factorises, and fails if the covariance is not
/// positive definite.
pub struct MvNormal {
  priv d : Int
  priv mean : FixedArray[Double]
  priv chol : FixedArray[Double]
}

///|
pub fn MvNormal::new(
  mean : ArrayView[Double],
  covariance : ArrayView[Double],
) -> MvNormal raise @dist.ParamError {
  let d = mean.length()
  if covariance.length() != d * d {
    raise @dist.InvalidParam(
      name="covariance",
      value=covariance.length().to_double(),
      reason="must be d * d, row-major",
    )
  }
  ...
}

///|
pub impl Dimension for MvNormal with dim(self) {
  self.d
}

///|
pub impl MvDensity for MvNormal with logpdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl MvSample for MvNormal with sample_into(self, rng, out) {
  ignore((self, rng, out))
  ...
}

///|
pub impl MvMoments for MvNormal with mean_into(self, out) {
  for i in 0..<self.d {
    out[i] = self.mean[i]
  }
  true
}

///|
pub impl MvMoments for MvNormal with covariance_into(self, out) {
  ignore((self, out))
  ...
}

///|
/// The Dirichlet distribution on the probability simplex.
pub struct Dirichlet {
  priv alpha : FixedArray[Double]
}

///|
pub fn Dirichlet::new(
  alpha : ArrayView[Double],
) -> Dirichlet raise @dist.ParamError {
  for a in alpha {
    if !(a > 0.0) {
      raise @dist.InvalidParam(name="alpha", value=a, reason="must be > 0")
    }
  }
  { alpha: FixedArray::from_iter(alpha.iter()) }
}

///|
pub impl Dimension for Dirichlet with dim(self) {
  self.alpha.length()
}

///|
pub impl MvDensity for Dirichlet with logpdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl MvSample for Dirichlet with sample_into(self, rng, out) {
  ignore((self, rng, out))
  ...
}

///|
/// Counts of n trials over d categories.
pub struct Multinomial {
  priv n : Int
  priv p : FixedArray[Double]
}

///|
pub fn Multinomial::new(
  n : Int,
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
pub impl Dimension for Multinomial with dim(self) {
  self.p.length()
}

///|
pub impl MvMass for Multinomial with logpmf(self, k) {
  ignore((self, k))
  ...
}

///|
/// A copula: a distribution on [0, 1]^d with uniform marginals. Combined
/// with any univariate marginals it builds a joint distribution whose
/// dependence is separate from its marginals (Sklar's theorem).
pub(open) trait Copula: Dimension {
  copula_cdf(Self, ArrayView[Double]) -> Double
  sample_uniform_into(Self, &@num.RandomSource, FixedArray[Double]) -> Unit
}

///|
/// The Gaussian copula with a given correlation matrix.
pub struct GaussianCopula {
  priv mvn : MvNormal
}

///|
pub fn GaussianCopula::new(
  correlation : ArrayView[Double],
) -> GaussianCopula raise @dist.ParamError {
  ignore(correlation)
  ...
}

///|
pub impl Dimension for GaussianCopula with dim(self) {
  self.mvn.d
}

///|
pub impl Copula for GaussianCopula with copula_cdf(self, u) {
  ignore((self, u))
  ...
}

///|
pub impl Copula for GaussianCopula with sample_uniform_into(self, rng, out) {
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
  { copula, marginals }
}

///|
pub impl[C : Copula] Dimension for Joint[C] with dim(self) {
  self.marginals.length()
}

///|
/// Sample uniforms from the copula, then map each through its marginal's
/// quantile function.
pub impl[C : Copula] MvSample for Joint[C] with sample_into(self, rng, out) {
  self.copula.sample_uniform_into(rng, out)
  for i, m in self.marginals {
    out[i] = @dist.Quantile::quantile(m, out[i])
  }
}
