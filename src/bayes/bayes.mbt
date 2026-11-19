///|
/// Computational methods beyond closed forms: MCMC and EM. Conjugate
/// updates live in `dist` (`Conjugate[P, D]`).

///|
/// Draws from a sampler, one row per draw, one column per parameter.
pub struct Chain {
  draws : @data.DataMatrix
  acceptance : Double
}

///|
/// Random-walk Metropolis–Hastings on an unnormalised log density.
pub fn metropolis(
  log_density : (ArrayView[Double]) -> Double,
  start : ArrayView[Double],
  step : ArrayView[Double],
  draws~ : Int32,
  burn_in~ : Int32,
  rng~ : &@num.RandomSource,
) -> Chain {
  ignore((log_density, start, step, draws, burn_in, rng))
  ...
}

///|
/// Gibbs sampling from full conditionals: `update(i, state, rng)` returns
/// a draw of coordinate i given the rest.
pub fn gibbs(
  update : (Int32, ArrayView[Double], &@num.RandomSource) -> Double,
  start : ArrayView[Double],
  draws~ : Int32,
  burn_in~ : Int32,
  rng~ : &@num.RandomSource,
) -> Chain {
  ignore((update, start, draws, burn_in, rng))
  ...
}

///|
/// Convergence diagnostics across chains for each parameter.
pub struct Diagnostics {
  /// Split R-hat (Vehtari et al. 2021).
  r_hat : FixedArray[Double]
  effective_size : FixedArray[Double]
}

///|
pub fn diagnose(chains : ArrayView[Chain]) -> Diagnostics {
  ignore(chains)
  ...
}

///|
/// Posterior summary of one parameter: mean, sd and a credible interval
/// (an `Interval` of kind `Credible`).
pub fn summarize(
  chain : Chain,
  param : Int32,
  level : Double,
) -> @infer.Estimate {
  ignore((chain, param, level))
  ...
}

///|
/// EM for a finite mixture of normals.
pub struct MixtureFit {
  weights : FixedArray[Double]
  means : FixedArray[Double]
  sds : FixedArray[Double]
  loglik : Double
  iterations : Int32
  converged : Bool
}

///|
pub fn em_normal_mixture(
  x : ArrayView[Double],
  components : Int32,
  rng~ : &@num.RandomSource,
  max_iterations~ : Int32 = 500,
) -> MixtureFit {
  ignore((x, components, rng, max_iterations))
  ...
}
