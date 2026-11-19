///|
/// Choices that R, SciPy, statsmodels and Julia resolve differently.
/// Each is a type the caller passes; none has a hidden default.

///|
/// How a two-sided p-value is formed when the null distribution is not
/// symmetric (binomial, hypergeometric, permutation nulls).
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum TwoSidedRule {
  /// Twice the smaller tail, capped at 1.
  Doubling
  /// The total probability of outcomes no more likely than the observed
  /// one (R's binom.test and fisher.test). Needs a mass function.
  MinLikelihood
}

///|
/// How the p-value was obtained. (The names end in `Null` so that they
/// cannot be taken for the cases of `Exactness`, which asks for a null.)
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum PValueKind {
  /// From the exact null distribution.
  ExactNull
  /// From a limiting distribution.
  AsymptoticNull
  /// From random draws of the null, with the (k + 1)/(N + 1) correction
  /// when `plus_one` is set, which keeps the p-value valid.
  MonteCarlo(draws~ : Int32, plus_one~ : Bool)
}

///|
/// When to use the exact null distribution. Libraries switch at different
/// sample sizes; here the threshold is stated.
#derive(Eq, Debug)
pub(all) enum Exactness {
  Exact
  Asymptotic
  /// Exact when the larger sample has at most this many observations and
  /// there are no ties.
  ExactUpTo(Int32)
}

///|
/// Pooled or separate variances. R's t.test defaults to separate (Welch),
/// SciPy's ttest_ind to pooled; this library has no default. Pooling
/// assumes equal variances, so that case carries the certificate: it
/// cannot be chosen without stating the assumption.
#derive(Eq, Debug)
pub(all) enum Variances {
  Pooled(Assumed[EqualVariance])
  Welch
}

///|
/// Zero differences in the Wilcoxon signed-rank test.
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum Zeros {
  /// Drop them (Wilcoxon's original).
  Drop
  /// Rank them, then drop their ranks (Pratt).
  Pratt
  /// Split their ranks between the signs.
  Split
}

///|
/// Tolerance for counting a resampled statistic as equal to the observed
/// one. Without it, rounding makes statistics that are equal in exact
/// arithmetic compare unequal and shifts randomisation p-values; SciPy
/// uses a relative tolerance for the same reason.
pub(all) struct TieTolerance {
  relative : Double
}
