///|
/// Resampling with user-defined statistics, after SciPy's bootstrap,
/// permutation_test and monte_carlo_test: the statistic is a function of
/// the sample, and the resampling scheme is a separate argument.

///|
/// A statistic of a sample of any shape.
pub type Statistic = (Sample) -> Double

///|
/// How a bootstrap interval is formed from the bootstrap distribution.
#derive(Eq, Debug)
pub(all) enum BootstrapInterval {
  Percentile
  Basic
  /// Bias-corrected and accelerated; needs a jackknife pass.
  BCa
}

///|
/// What is exchanged under the null in a permutation test.
#derive(Eq, Debug)
pub(all) enum PermutationScheme {
  /// Observations are exchanged between independent samples.
  Independent
  /// The pairing between paired observations is broken.
  Pairings
  /// Within each pair, the two values are swapped (sign flips).
  WithinPairs
}

///|
pub struct BootstrapResult {
  estimate : Double
  interval : Interval
  se : Double
  /// The bootstrap distribution itself, as an ordinary distribution.
  distribution : @dist.Empirical
  reps : Int
}

///|
/// Resamples each independent sample with replacement (pairs together for
/// `Paired` and `Bivariate`-like shapes).
pub fn bootstrap(
  sample : Sample,
  statistic : Statistic,
  reps~ : Int,
  level~ : Double,
  interval~ : BootstrapInterval,
  rng~ : &@num.RandomSource,
) -> BootstrapResult raise QuestionError {
  ignore((sample, statistic, reps, level, interval, rng))
  ...
}

///|
/// The permutation null of a statistic. `reps = None` enumerates every
/// permutation (exact); otherwise draws that many, and the p-values carry
/// the (k + 1)/(N + 1) correction.
pub fn permutation_test(
  sample : Sample,
  statistic : Statistic,
  scheme~ : PermutationScheme,
  reps~ : Int?,
  ties~ : TieTolerance,
  rng~ : &@num.RandomSource,
) -> NullModel raise QuestionError {
  ignore((sample, statistic, scheme, reps, ties, rng))
  ...
}

///|
/// A Monte Carlo test: the null is a generator of samples of the same
/// shape, and the statistic is computed on each draw.
pub fn monte_carlo_test(
  observed : Sample,
  null_sampler : (&@num.RandomSource) -> Sample,
  statistic : Statistic,
  reps~ : Int,
  rng~ : &@num.RandomSource,
) -> NullModel raise QuestionError {
  ignore((observed, null_sampler, statistic, reps, rng))
  ...
}

///|
/// Leave-one-out estimates of bias and standard error.
pub fn jackknife(sample : Sample, statistic : Statistic) -> Estimate raise QuestionError {
  ignore((sample, statistic))
  ...
}
