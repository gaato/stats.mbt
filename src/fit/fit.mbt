///|
/// Fitting families to data. Every fit returns the distribution together
/// with what is needed to judge it: the log-likelihood, the number of free
/// parameters and observations, whether the optimiser converged, and the
/// asymptotic covariance when it is available.

///|
#derive(Debug)
pub struct MleResult[D] {
  dist : D
  loglik : Double
  /// Free parameters of the family.
  k : Int32
  /// Observations.
  n : Int64
  /// k × k, on the family's own parameter scale; `None` when the
  /// information matrix is singular or was not computed.
  covariance : @data.DataMatrix?
  converged : Bool
  iterations : Int32
  family : String
}

///|
/// Akaike's information criterion, 2k - 2 ln L.
pub fn[D] MleResult::aic(self : MleResult[D]) -> Double {
  2.0 * self.k.to_double() - 2.0 * self.loglik
}

///|
/// The small-sample correction; `None` when n <= k + 1 makes it undefined.
pub fn[D] MleResult::aicc(self : MleResult[D]) -> Double? {
  let n = self.n.to_double()
  let k = self.k.to_double()
  if n <= k + 1.0 {
    None
  } else {
    Some(self.aic() + 2.0 * k * (k + 1.0) / (n - k - 1.0))
  }
}

///|
/// Schwarz's Bayesian information criterion, k ln n - 2 ln L.
pub fn[D] MleResult::bic(self : MleResult[D]) -> Double {
  self.k.to_double() * @math.ln(self.n.to_double()) - 2.0 * self.loglik
}

///|
/// Maximum likelihood for an exponential family: the estimate is the
/// member whose expected sufficient statistic equals the sample average.
/// No optimiser is involved, so the result always reports convergence.
pub fn[D : @dist.Density] fit_exp_family(
  ef : @dist.ExpFamily[D],
  data : ArrayView[Double],
) -> MleResult[D] raise @dist.FitError {
  let k = ef.family.dim()
  if data.length() <= k {
    raise @dist.TooFewPoints(needed=k + 1, got=data.length())
  }
  let mean = FixedArray::make(ef.natural_dim, 0.0)
  for x in data {
    let t = (ef.sufficient)(x)
    for j in 0..<ef.natural_dim {
      mean[j] = mean[j] + t[j]
    }
  }
  let n = data.length().to_int64()
  for j in 0..<ef.natural_dim {
    mean[j] = mean[j] / n.to_double()
  }
  let d = (ef.from_mean_sufficient)(mean[:])
  {
    dist: d,
    loglik: @dist.loglik(d, data),
    k,
    n,
    covariance: None,
    converged: true,
    iterations: 0,
    family: ef.family.name,
  }
}

///|
/// General maximum likelihood: optimise on the unconstrained scale from
/// `start`, using the family's closed-form score.
pub fn[D : @dist.Density + @dist.Score] fit_mle(
  family : @dist.Family[D],
  data : ArrayView[Double],
  start~ : ArrayView[Double],
  max_iterations~ : Int32 = 200,
) -> MleResult[D] raise @dist.FitError {
  ignore((family, data, start, max_iterations))
  ...
}

///|
/// Pick among fitted candidates by AIC. The candidates may come from
/// different families, so they are passed as already-fitted summaries.
#derive(Eq, Debug)
pub(all) struct Candidate {
  family : String
  aic : Double
  bic : Double
  k : Int32
}

///|
pub fn[D] MleResult::candidate(self : MleResult[D]) -> Candidate {
  { family: self.family, aic: self.aic(), bic: self.bic(), k: self.k, }
}

///|
/// The candidate with the smallest AIC and the Akaike weight of each:
/// exp(-Δᵢ/2) normalised, where Δᵢ is the AIC difference from the best.
pub fn akaike_weights(cands : ArrayView[Candidate]) -> FixedArray[Double] {
  let w = FixedArray::make(cands.length(), 0.0)
  if cands.length() == 0 {
    return w
  }
  let mut best = cands[0].aic
  for c in cands {
    if c.aic < best {
      best = c.aic
    }
  }
  let mut total = 0.0
  for i, c in cands {
    w[i] = @math.exp(-0.5 * (c.aic - best))
    total = total + w[i]
  }
  for i in 0..<w.length() {
    w[i] = w[i] / total
  }
  w
}
