///|
/// Results of inference. Tests return a value that keeps the statistic,
/// the reference distribution and the settings together, so a p-value
/// never travels alone.

///|
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
/// assumption therefore always appears in the caller's code.
pub struct Assumed[A] {
  priv _phantom : Unit
}

///|
pub fn[A] Assumed::postulate() -> Assumed[A] {
  { _phantom: () }
}

///|
pub(all) enum Normality {}

///|
pub(all) enum EqualVariance {}

///|
/// What a test returns. Fields follow the union of R's htest, SciPy's
/// result objects and HypothesisTests.jl; anything a test cannot provide is
/// `None` rather than a placeholder. The null model is kept, so the
/// p-value for another alternative needs no recomputation.
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
/// A confidence interval carries its level and the method that produced
/// it; there is no bare (lo, hi) pair in the public API.
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
pub struct Estimate {
  value : Double
  se : Double?
  interval : Interval?
  n : Int64
  procedure : String
}

///|
/// The README-first check: a t test written the way a user would write it.
/// Type-checked only; the bodies it calls are still `...`.
fn _usage_sketch() -> Unit {
  let a = [5.1, 4.9, 5.6, 5.8, 6.0]
  let b = [6.2, 6.4, 5.9, 6.8, 7.1]
  // The assumption is visible at the call site.
  let normality : Assumed[Normality] = Assumed::postulate()
  let r = try
    Some(
      t_test_two(
        a[:],
        b[:],
        variances=Welch,
        alternative=TwoSided,
        level=0.95,
        normality~,
      ),
    )
  catch {
    _ => None
  }
  ignore(r)
  // Summaries compose: the two halves are summarized separately.
  let empty : @describe.MeanVar = @describe.Accumulator::empty()
  let left = a.fold(init=empty, (acc, x) => acc.push(x))
  let right = b.fold(init=empty, (acc, x) => acc.push(x))
  let both = left.merge(right)
  ignore(both.variance(ddof=1))
  // Reference distributions are ordinary distributions.
  let t = try Some(@dist.StudentT::new(8.0)) catch { _ => None }
  ignore(t)
}
