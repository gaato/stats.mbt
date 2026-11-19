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
#non_exhaustive
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
/// (Named in the plural: `@dist.Sample` is the trait for drawing from a
/// distribution.)
pub(all) enum Samples {
  One(ArrayView[Double])
  Two(ArrayView[Double], ArrayView[Double])
  /// Two measurements on the same units, or (x, y) pairs. A `Bivariate`
  /// has checked that the two have the same length.
  Paired(@data.Bivariate)
  /// k independent groups.
  Groups(Array[ArrayView[Double]])
  /// k related measurements on each of n blocks.
  Blocks(@data.Blocks)
}

///|
pub struct Question {
  estimand : Estimand
  sample : Samples
}

///|
/// Why a question could not be answered. The cases separate what the
/// caller can fix (shapes, lengths, arguments) from what the data make
/// impossible and from what this library does not do.
#non_exhaustive
#derive(Debug)
pub(all) suberror QuestionError {
  /// The estimand does not fit the data shape (a mean difference from one
  /// sample), so the mismatch is an error at the point of asking, not a
  /// NaN later.
  ShapeMismatch(estimand~ : Estimand)
  /// The permutation scheme does not fit the data shape (pairings of two
  /// independent samples).
  SchemeMismatch(scheme~ : PermutationScheme)
  /// Two arrays that must have the same length do not.
  LengthMismatch(expected~ : Int32, got~ : Int32)
  TooFewPoints(needed~ : Int32, got~ : Int32)
  /// An argument outside what the procedure accepts (a level outside
  /// (0, 1), fewer than two replicates).
  InvalidArgument(name~ : String, reason~ : String)
  /// The data make the statistic undefined (zero variance, empty cell).
  Degenerate(reason~ : String)
  /// The request has no answer by the method asked for, however it is
  /// implemented (an interval from a permutation; more arrangements than
  /// can be enumerated).
  Unsupported(reason~ : String)
  /// A procedure or path this library does not compute yet. It is named,
  /// never approximated silently.
  NotImplemented(what~ : String)
}

///|
/// Checks the estimand against the shape: one number must answer it. A
/// mean of several groups is not one number, so groups and blocks are
/// taken by the tests that compare them, not by a question.
pub fn Question::new(
  estimand : Estimand,
  sample : Samples,
) -> Question raise QuestionError {
  let ok = match (estimand, sample) {
    (Mean | Median | Quantile(_) | Variance | Proportion, One(_)) => true
    (MeanDifference | MedianDifference | MeanRatio, Two(_, _) | Paired(_)) =>
      true
    (Correlation, Paired(_)) => true
    _ => false
  }
  if !ok {
    raise ShapeMismatch(estimand~)
  }
  if estimand is Quantile(p) && !(p >= 0.0 && p <= 1.0) {
    raise InvalidArgument(name="p", reason="a quantile is asked at p in [0, 1]")
  }
  { estimand, sample, }
}

///|
/// How many resamples a randomisation method uses.
#derive(Eq, Debug)
pub(all) enum Resamples {
  /// Every arrangement, enumerated: the exact randomisation distribution.
  All
  /// This many random arrangements; p-values then carry the
  /// (k + 1)/(N + 1) correction.
  Random(Int32)
}

///|
/// How uncertainty is assessed. Resampling methods take their own
/// replicate count; the random source is passed separately.
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum Method {
  /// The textbook sampling distribution, under the assumptions it needs.
  Analytic
  Bootstrap(reps~ : Int32)
  Permutation(Resamples)
}

///|
/// The null hypothesis.
#non_exhaustive
#derive(Eq, Debug)
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
  check_level(level)
  match how {
    Analytic => analytic_estimate(q, level)
    Bootstrap(reps~) => {
      let stat = estimand_statistic(q.estimand)
      let b = bootstrap(
        q.sample,
        stat,
        reps~,
        level~,
        interval=Percentile,
        rng~,
      )
      {
        value: b.estimate,
        se: Some(b.se),
        interval: Some(b.interval),
        n: sample_size(q.sample),
        procedure: "bootstrap (percentile)",
      }
    }
    Permutation(_) =>
      raise Unsupported(
        reason="permutation gives tests, not intervals; use null_model",
      )
  }
}

///|
fn sample_size(s : Samples) -> Int64 {
  match s {
    One(x) => x.length().to_int64()
    Two(x, y) => (x.length() + y.length()).to_int64()
    Paired(xy) => xy.length().to_int64()
    Groups(gs) => {
      let mut n = 0
      for g in gs {
        n = n + g.length()
      }
      n.to_int64()
    }
    Blocks(b) => (b.blocks * b.treatments).to_int64()
  }
}

///|
fn mean_of(x : ArrayView[Double]) -> Double {
  summarize(x).1
}

///|
fn median_of(x : ArrayView[Double]) -> Double {
  @describe.quantile(x, 0.5, Type7) catch {
    _ => 0.0 / 0.0
  }
}

///|
fn differences(xy : @data.Bivariate) -> FixedArray[Double] {
  FixedArray::makei(xy.length(), i => xy.x[i] - xy.y[i])
}

///|
/// The statistic an estimand stands for, on the sample shapes it accepts.
fn estimand_statistic(e : Estimand) -> Statistic {
  fn(s : Samples) -> Double {
    match (e, s) {
      (Mean | Proportion, One(x)) => mean_of(x)
      (Median, One(x)) => median_of(x)
      (Quantile(p), One(x)) =>
        @describe.quantile(x, p, Type7) catch {
          _ => 0.0 / 0.0
        }
      (Variance, One(x)) => summarize(x).2
      (MeanDifference, Two(x, y)) => mean_of(x) - mean_of(y)
      (MeanDifference, Paired(xy)) => mean_of(differences(xy)[:])
      (MedianDifference, Two(x, y)) => median_of(x) - median_of(y)
      (MedianDifference, Paired(xy)) => median_of(differences(xy)[:])
      (MeanRatio, Two(x, y)) => mean_of(x) / mean_of(y)
      (MeanRatio, Paired(xy)) => mean_of(xy.x[:]) / mean_of(xy.y[:])
      (Correlation, Paired(xy)) => {
        let mut c = @describe.Comoment::empty()
        for i in 0..<xy.length() {
          c = c.push(xy.x[i], xy.y[i])
        }
        match c.correlation() {
          Some(r) => r
          None => 0.0 / 0.0
        }
      }
      // `Question::new` admits no other pair
      _ => 0.0 / 0.0
    }
  }
}

///|
/// The textbook intervals: t for means and mean differences, χ² for a
/// variance, Clopper–Pearson for a proportion, Fisher z for a correlation.
fn analytic_estimate(
  q : Question,
  level : Double,
) -> Estimate raise QuestionError {
  let normality : Assumed[Normality] = Assumed::postulate()
  let from_test = (r : TestResult, n : Int64, procedure : String) => {
    let se = match (r.estimate, r.null_value) {
      (Some(e), Some(v)) if r.statistic != 0.0 =>
        Some(((e - v) / r.statistic).abs())
      _ => None
    }
    { value: r.estimate.unwrap(), se, interval: r.interval, n, procedure, }
  }
  match (q.estimand, q.sample) {
    (Mean, One(x)) => {
      let r = t_test_one(x, mu=0.0, alternative=TwoSided, level~, normality~)
      from_test(r, x.length().to_int64(), "t interval")
    }
    (MeanDifference, Two(x, y)) => {
      let r = t_test_two(
        x,
        y,
        variances=Welch,
        alternative=TwoSided,
        level~,
        normality~,
      )
      from_test(r, (x.length() + y.length()).to_int64(), "Welch t interval")
    }
    (MeanDifference, Paired(xy)) => {
      let r = t_test_paired(xy, alternative=TwoSided, level~, normality~)
      from_test(r, xy.length().to_int64(), "paired t interval")
    }
    (Variance, One(x)) => {
      let (n, _, v) = summarize(x)
      if n < 2 {
        raise TooFewPoints(needed=2, got=n)
      }
      let df = (n - 1).to_double()
      let chi = @dist.ChiSquared::new(df) catch { _ => panic() }
      let alpha = 1.0 - level
      {
        value: v,
        se: Some(v * (2.0 / df).sqrt()),
        interval: Some({
          lo: df * v / chi.isf(0.5 * alpha),
          hi: df * v / chi.quantile(0.5 * alpha),
          level,
          kind: Confidence,
          procedure: "chi-square interval",
        }),
        n: n.to_int64(),
        procedure: "chi-square interval",
      }
    }
    (Proportion, One(x)) => {
      let mut k = 0L
      for v in x {
        if v == 1.0 {
          k = k + 1L
        } else if v != 0.0 {
          raise Degenerate(reason="a proportion needs 0/1 data")
        }
      }
      let n = x.length().to_int64()
      if n == 0L {
        raise TooFewPoints(needed=1, got=0)
      }
      let p = k.to_double() / n.to_double()
      {
        value: p,
        se: Some((p * (1.0 - p) / n.to_double()).sqrt()),
        interval: Some(clopper_pearson(k, n, level, TwoSided)),
        n,
        procedure: "Clopper-Pearson",
      }
    }
    (Correlation, Paired(xy)) => {
      let r = correlation_test(xy, Pearson, alternative=TwoSided, level~)
      {
        value: r.estimate.unwrap(),
        se: None,
        interval: r.interval,
        n: xy.length().to_int64(),
        procedure: "Fisher z interval",
      }
    }
    _ =>
      raise Unsupported(
        reason="no analytic interval for this estimand; use Bootstrap",
      )
  }
}

///|
/// The null distribution of the question's statistic under the hypothesis.
pub fn null_model(
  q : Question,
  h : Hypothesis,
  how : Method,
  rng~ : &@num.RandomSource,
) -> NullModel raise QuestionError {
  let stat = estimand_statistic(q.estimand)
  match (how, h, q.sample) {
    (Permutation(reps), Independence, Two(_, _)) =>
      permutation_test(
        q.sample,
        stat,
        scheme=Independent,
        reps~,
        ties={ relative: 1.0e-12, },
        rng~,
      )
    (Permutation(reps), Independence, Paired(_)) =>
      permutation_test(
        q.sample,
        stat,
        scheme=if q.estimand == Correlation { Pairings } else { WithinPairs },
        reps~,
        ties={ relative: 1.0e-12, },
        rng~,
      )
    _ =>
      raise NotImplemented(
        what="null_model beyond permutation nulls of independence; the analytic tests are separate functions",
      )
  }
}
