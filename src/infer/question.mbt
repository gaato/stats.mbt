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
  TooFewPoints(needed~ : Int32, got~ : Int32)
  /// A procedure that relies on an assumption was called without its
  /// certificate (a pooled t test without `Assumed[EqualVariance]`).
  MissingAssumption(name~ : String)
  /// The data make the statistic undefined (zero variance, empty cell).
  Degenerate(reason~ : String)
  /// A combination this library does not compute (yet), named rather than
  /// approximated silently.
  Unsupported(reason~ : String)
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
  { estimand, sample, }
}

///|
/// How uncertainty is assessed. Resampling methods take their own
/// replicate count; the random source is passed separately.
pub(all) enum Method {
  /// The textbook sampling distribution, under the assumptions it needs.
  Analytic
  Bootstrap(reps~ : Int32)
  Permutation(reps~ : Int32)
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
    Permutation(..) =>
      raise Unsupported(
        reason="permutation gives tests, not intervals; use null_model",
      )
  }
}

///|
fn sample_size(s : Sample) -> Int64 {
  match s {
    One(x) => x.length().to_int64()
    Two(x, y) | Paired(x, y) => (x.length() + y.length()).to_int64()
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
/// The statistic an estimand stands for, on the sample shapes it accepts.
fn estimand_statistic(e : Estimand) -> Statistic {
  fn(s : Sample) -> Double {
    match (e, s) {
      (Mean | Proportion, One(x)) => mean_of(x)
      (Median, One(x)) => median_of(x)
      (Quantile(p), One(x)) =>
        @describe.quantile(x, p, Type7) catch {
          _ => 0.0 / 0.0
        }
      (Variance, One(x)) => summarize(x).2
      (MeanDifference, Two(x, y)) => mean_of(x) - mean_of(y)
      (MeanDifference, Paired(x, y)) => {
        let d = FixedArray::makei(x.length(), i => x[i] - y[i])
        mean_of(d[:])
      }
      (MedianDifference, Two(x, y)) => median_of(x) - median_of(y)
      (MedianDifference, Paired(x, y)) => {
        let d = FixedArray::makei(x.length(), i => x[i] - y[i])
        median_of(d[:])
      }
      (MeanRatio, Two(x, y) | Paired(x, y)) => mean_of(x) / mean_of(y)
      (Correlation, Paired(x, y)) => {
        let mut c = @describe.Comoment::empty()
        for i in 0..<x.length() {
          c = c.push(x[i], y[i])
        }
        match c.correlation() {
          Some(r) => r
          None => 0.0 / 0.0
        }
      }
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
    (MeanDifference, Paired(x, y)) => {
      let r = t_test_paired(x, y, alternative=TwoSided, level~, normality~)
      from_test(r, x.length().to_int64(), "paired t interval")
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
          lo: df * v / @dist.Quantile::isf(chi, 0.5 * alpha),
          hi: df * v / @dist.Quantile::quantile(chi, 0.5 * alpha),
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
    (Correlation, Paired(x, y)) => {
      let xy = @data.Bivariate::new(x, y) catch {
        _ => raise ShapeMismatch(estimand=Correlation)
      }
      let r = correlation_test(xy, Pearson, alternative=TwoSided, level~)
      {
        value: r.estimate.unwrap(),
        se: None,
        interval: r.interval,
        n: x.length().to_int64(),
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
  let stat = estimand_statistic(q.estimand)
  match (how, h, q.sample) {
    (Permutation(reps~), Independence, Two(_, _)) =>
      permutation_test(
        q.sample,
        stat,
        scheme=Independent,
        reps=Some(reps),
        ties={ relative: 1.0e-12, },
        rng~,
      )
    (Permutation(reps~), Independence, Paired(_, _)) =>
      permutation_test(
        q.sample,
        stat,
        scheme=if q.estimand == Correlation { Pairings } else { WithinPairs },
        reps=Some(reps),
        ties={ relative: 1.0e-12, },
        rng~,
      )
    _ =>
      raise Unsupported(
        reason="null_model supports permutation nulls of independence for now; the analytic tests are separate functions",
      )
  }
}

///|
/// The p-value from a null model, the same code for theory-based and
/// randomisation tests. Each tail includes the observed value itself
/// (P(T >= t), P(T <= t)), which matters for discrete and permutation nulls,
/// and neither tail is computed as 1 minus the other.
pub fn NullModel::p_value(
  self : NullModel,
  alternative : Alternative,
) -> Double {
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
        // needs the null's mass function, which `Univariate` does not carry;
        // tests that use this rule compute it from their discrete null
        MinLikelihood => 0.0 / 0.0
      }
  }
}
