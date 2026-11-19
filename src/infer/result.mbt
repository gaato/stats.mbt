///|
/// Results of inference. Tests return a value that keeps the statistic,
/// the reference distribution and the settings together, so a p-value
/// never travels alone.

///|
/// The alternative hypothesis: the parameter differs from, is below, or is
/// above its null value.
#derive(Eq, Debug)
pub(all) enum Alternative {
  TwoSided
  Less
  Greater
}

///|
/// Assumptions as certificates, in the style of Luna-Flow's `Law[Th, A]`:
/// a test that relies on normality asks for an `Assumed[Normality]`, and the
/// only ways to get one are a diagnostic or an explicit `postulate()`. The
/// assumption therefore always appears in the caller's code. A choice that
/// brings an assumption with it carries the certificate in the choice
/// itself (`Variances::Pooled`), so it cannot be made without one.
#derive(Eq, Debug)
pub struct Assumed[A] {
  priv _phantom : Unit
}

///|
/// The caller's statement that the assumption holds, made explicitly.
pub fn[A] Assumed::postulate() -> Assumed[A] {
  { _phantom: (), }
}

///|
/// The assumption that the data come from a normal distribution.
#derive(Eq, Debug)
pub(all) enum Normality {}

///|
/// The assumption that the groups share one variance.
#derive(Eq, Debug)
pub(all) enum EqualVariance {}

///|
/// What a test returns. Fields follow the union of R's htest, SciPy's
/// result objects and HypothesisTests.jl; anything a test cannot provide is
/// `None` rather than a placeholder. The null model is kept, so the
/// p-value for another alternative needs no recomputation.
#derive(Debug)
pub struct TestResult {
  statistic : Double
  p_value : Double
  /// Log of the p-value, for values below the smallest positive Double.
  log_p_value : Double
  df : Double?
  /// Second degrees of freedom (F tests).
  df2 : Double?
  alternative : Alternative
  /// The estimate the test is about (a mean difference, an odds ratio).
  estimate : Double?
  /// Its value under the null hypothesis.
  null_value : Double?
  interval : Interval?
  effect : Effect?
  p_kind : PValueKind
  null_model : NullModel?
  procedure : String
}

///|
/// One comparison among several groups: which two (by their position in
/// the input, `first < second`) and the test between them.
#derive(Debug)
pub struct Comparison {
  first : Int32
  second : Int32
  result : TestResult
}

///|
/// A confidence interval carries its level and the method that produced
/// it; there is no bare (lo, hi) pair in the public API.
#derive(Eq, Debug)
pub struct Interval {
  lo : Double
  hi : Double
  level : Double
  kind : IntervalKind
  procedure : String
}

///|
/// Intervals with different meanings are different values of this, never
/// the same struct with a comment.
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum IntervalKind {
  /// Covers the true parameter with the stated frequency.
  Confidence
  /// Holds the stated posterior probability.
  Credible
  /// Covers a future observation.
  Prediction
  /// Covers the stated proportion of the population.
  Tolerance
}

///|
/// Point estimate with optional uncertainty. `se` and `interval` are
/// `None` when the method gives none, rather than a misleading 0 or NaN.
#derive(Eq, Debug)
pub struct Estimate {
  value : Double
  se : Double?
  interval : Interval?
  n : Int64
  procedure : String
}

///|
/// The README-first check: a t test written the way a user would write it.
/// Type-checked only.
fn _usage_sketch() -> Unit {
  let a = [5.1, 4.9, 5.6, 5.8, 6.0]
  let b = [6.2, 6.4, 5.9, 6.8, 7.1]
  // The assumption is visible at the call site.
  let normality : Assumed[Normality] = Assumed::postulate()
  let r = Some(
    t_test_two(
      a[:],
      b[:],
      variances=Welch,
      alternative=TwoSided,
      level=0.95,
      normality~,
    ),
  ) catch {
    _ => None
  }
  ignore(r)
  // Summaries compose: the two halves are summarized separately.
  let left : @describe.MeanVar = @describe.accumulate(a)
  let right : @describe.MeanVar = @describe.accumulate(b)
  let both = left.merge(right)
  ignore(both.variance(ddof=1))
  // Reference distributions are ordinary distributions.
  let t = Some(@dist.StudentT::new(8.0)) catch { _ => None }
  ignore(t)
}
