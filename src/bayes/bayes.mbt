///|
/// Computational methods beyond closed forms: MCMC and EM. Conjugate
/// updates live in `dist` (`Conjugate[P, D]`). Failures are
/// `@regress.ModelError`.

///|
/// Draws from a sampler, one row per draw, one column per parameter.
#derive(Eq, Debug)
pub struct Chain {
  draws : @data.DataMatrix
  /// The share of proposals accepted; `None` for draws not made by a
  /// Metropolis step (Gibbs updates, draws from elsewhere).
  acceptance : Double?
}

///|
/// A chain from draws made elsewhere (another sampler, a file), so that
/// `diagnose` and `summarize` work on them.
pub fn Chain::from_draws(draws : @data.DataMatrix) -> Chain {
  { draws, acceptance: None, }
}

///|
/// Random-walk Metropolis–Hastings on an unnormalised log density.
pub fn metropolis(
  log_density : (ArrayView[Double]) -> Double,
  start~ : ArrayView[Double],
  step~ : ArrayView[Double],
  draws~ : Int32,
  burn_in~ : Int32,
  rng~ : &@num.RandomSource,
) -> Chain raise @regress.ModelError {
  if start.length() != step.length() {
    raise LengthMismatch(expected=start.length(), got=step.length())
  }
  ignore((log_density, draws, burn_in, rng))
  ...
}

///|
/// Gibbs sampling from full conditionals: `update(i, state, rng)` returns
/// a draw of coordinate i given the rest.
pub fn gibbs(
  update : (Int32, ArrayView[Double], &@num.RandomSource) -> Double,
  start~ : ArrayView[Double],
  draws~ : Int32,
  burn_in~ : Int32,
  rng~ : &@num.RandomSource,
) -> Chain raise @regress.ModelError {
  ignore((update, start, draws, burn_in, rng))
  ...
}

///|
/// Convergence diagnostics across chains for each parameter.
#derive(Eq, Debug)
pub struct Diagnostics {
  /// Split R-hat (Vehtari et al. 2021).
  r_hat : FixedArray[Double]
  effective_size : FixedArray[Double]
}

///|
/// Every chain must have the same parameters (columns).
pub fn diagnose(
  chains : ArrayView[Chain],
) -> Diagnostics raise @regress.ModelError {
  ignore(chains)
  ...
}

///|
/// Posterior summary of one parameter: mean, sd and a credible interval
/// (an `Interval` of kind `Credible`).
pub fn summarize(
  chain : Chain,
  parameter~ : Int32,
  level~ : Double,
) -> @infer.Estimate raise @regress.ModelError {
  ignore((chain, parameter, level))
  ...
}

///|
/// EM for a finite mixture of normals.
#derive(Eq, Debug)
pub struct MixtureFit {
  weights : FixedArray[Double]
  means : FixedArray[Double]
  sds : FixedArray[Double]
  loglik : Double
  iterations : Int32
  converged : Bool
}

///|
/// Fits a mixture of `components` normals to `x` by EM, from starting
/// points drawn with `rng`, for at most `max_iterations` steps.
pub fn em_normal_mixture(
  x : ArrayView[Double],
  components~ : Int32,
  rng~ : &@num.RandomSource,
  max_iterations~ : Int32 = 500,
) -> MixtureFit raise @regress.ModelError {
  ignore((x, components, rng, max_iterations))
  ...
}
