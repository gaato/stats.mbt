///|
/// The question-first layer, after R's infer package: say what is to be
/// estimated, from which data, by which method; the answer comes back with
/// its uncertainty. The algorithm is a choice made after the question, so
/// swapping a bootstrap for an analytic standard error does not change the
/// question.
///
/// The unifying observation: a test needs an observed statistic and a null
/// distribution, and a null distribution is an ordinary distribution. A
/// theory-based test uses a parametric one (Student t); a randomisation
/// test builds an `Empirical` one from permutations or bootstrap draws.
/// Either way the p-value is the null distribution's upper tail at the
/// observed statistic.

///|
/// What is to be estimated.
#derive(Eq, Debug)
pub(all) enum Estimand {
  Mean
  Median
  Quantile(Double)
  Variance
  Proportion
  MeanDifference
  MedianDifference
  MeanRatio
  Correlation
}

///|
/// The shape of the data, which fixes which estimands make sense.
/// Fixed now, because adding a case later would break callers'
/// exhaustive matches; the cases cover the shapes found in R, SciPy,
/// statsmodels and the mooncakes survey that resampling can act on.
pub(all) enum Sample {
  One(ArrayView[Double])
  Two(ArrayView[Double], ArrayView[Double])
  /// Two measurements on the same units, or (x, y) pairs.
  Paired(ArrayView[Double], ArrayView[Double])
  /// k independent groups.
  Groups(Array[ArrayView[Double]])
  /// k related measurements on each of n blocks.
  Blocks(@data.Blocks)
}

///|
pub struct Question {
  estimand : Estimand
  sample : Sample
}

///|
/// Raised when the estimand does not fit the data shape (a mean
/// difference from one sample), so the mismatch is an error at the point
/// of asking, not a NaN later.
#derive(Debug)
pub(all) suberror QuestionError {
  ShapeMismatch(estimand~ : Estimand)
  TooFewPoints(needed~ : Int, got~ : Int)
}

///|
pub fn Question::new(
  estimand : Estimand,
  sample : Sample,
) -> Question raise QuestionError {
  let ok = match (estimand, sample) {
    (Mean | Median | Quantile(_) | Variance | Proportion, One(_)) => true
    (MeanDifference | MedianDifference | MeanRatio, Two(_, _) | Paired(_, _)) =>
      true
    (Correlation, Paired(_, _)) => true
    (Mean | Median | Variance, Groups(_) | Blocks(_)) => true
    _ => false
  }
  if !ok {
    raise ShapeMismatch(estimand~)
  }
  { estimand, sample }
}

///|
/// How uncertainty is assessed. Resampling methods take their own
/// replicate count; the random source is passed separately.
pub(all) enum Method {
  /// The textbook sampling distribution, under the assumptions it needs.
  Analytic
  Bootstrap(reps~ : Int)
  Permutation(reps~ : Int)
}

///|
/// The null hypothesis.
pub(all) enum Hypothesis {
  /// The estimand equals a value.
  Point(Double)
  /// The two groups (or the paired variables) are independent.
  Independence
}

///|
/// Point estimate and interval for the question.
pub fn estimate(
  q : Question,
  how : Method,
  level~ : Double,
  rng~ : &@num.RandomSource,
) -> Estimate raise QuestionError {
  ignore((q, how, level, rng))
  ...
}

///|
/// Observed statistic and null distribution, before any p-value. Keeping
/// the null distribution lets a caller look at it, plot it, or compute a
/// p-value for another alternative without recomputing.
pub(all) struct NullModel {
  observed : Double
  null_dist : &@dist.Univariate
  two_sided : TwoSidedRule
  p_kind : PValueKind
  procedure : String
}

///|
pub fn null_model(
  q : Question,
  h : Hypothesis,
  how : Method,
  rng~ : &@num.RandomSource,
) -> NullModel raise QuestionError {
  ignore((q, h, how, rng))
  ...
}

///|
/// The p-value from a null model, the same code for theory-based and
/// randomisation tests. Each tail includes the observed value itself
/// (P(T >= t), P(T <= t)), which matters for discrete and permutation nulls,
/// and neither tail is computed as 1 minus the other.
pub fn NullModel::p_value(self : NullModel, alternative : Alternative) -> Double {
  let upper = @dist.Cdf::sf_at_least(self.null_dist, self.observed)
  let lower = @dist.Cdf::cdf(self.null_dist, self.observed)
  match alternative {
    Greater => upper
    Less => lower
    TwoSided =>
      match self.two_sided {
        Doubling => {
          let p = 2.0 * (if upper < lower { upper } else { lower })
          if p > 1.0 {
            1.0
          } else {
            p
          }
        }
        // needs the null's mass function, which `Univariate` does not carry
        MinLikelihood => ...
      }
  }
}
