///|
/// Shared pieces of the tests: two-pass summaries, critical values and
/// result assembly.

///|
/// (n, mean, variance with divisor n - 1), two-pass with compensated sums.
fn summarize(x : ArrayView[Double]) -> (Int32, Double, Double) {
  let n = x.length()
  let mut s = @num.CompensatedSum::zero()
  for v in x {
    s = s.add(v)
  }
  let m = s.value() / n.to_double()
  let mut ss = @num.CompensatedSum::zero()
  for v in x {
    ss = ss.add((v - m) * (v - m))
  }
  let v = if n > 1 { ss.value() / (n - 1).to_double() } else { 0.0 / 0.0 }
  (n, m, v)
}

///|
fn check_level(level : Double) -> Unit raise QuestionError {
  if !(level > 0.0 && level < 1.0) {
    raise InvalidArgument(name="level", reason="must be in (0, 1)")
  }
}

///|
/// The interval est ± crit·se in the direction of the alternative, with
/// crit the upper quantile of `pivot` at the tail probability required.
fn pivot_interval(
  pivot : &@dist.Univariate,
  est : Double,
  se : Double,
  level : Double,
  alternative : Alternative,
  procedure : String,
) -> Interval {
  let alpha = 1.0 - level
  match alternative {
    TwoSided => {
      let c = @dist.Quantile::isf(pivot, 0.5 * alpha)
      {
        lo: est - c * se,
        hi: est + c * se,
        level,
        kind: Confidence,
        procedure,
      }
    }
    Greater => {
      let c = @dist.Quantile::isf(pivot, alpha)
      { lo: est - c * se, hi: 1.0 / 0.0, level, kind: Confidence, procedure, }
    }
    Less => {
      let c = @dist.Quantile::isf(pivot, alpha)
      { lo: -1.0 / 0.0, hi: est + c * se, level, kind: Confidence, procedure, }
    }
  }
}

///|
/// The log of a tail probability of a continuous null, through `logsf` so
/// that p-values below 1e-308 are still reported on the log scale.
fn log_tail(
  null_dist : &@dist.Univariate,
  observed : Double,
  alternative : Alternative,
) -> Double {
  match alternative {
    Greater => @dist.Cdf::logsf(null_dist, observed)
    Less => @dist.Cdf::logcdf(null_dist, observed)
    TwoSided => {
      let lu = @dist.Cdf::logsf(null_dist, observed)
      let ll = @dist.Cdf::logcdf(null_dist, observed)
      let l = @math.ln(2.0) + (if lu < ll { lu } else { ll })
      if l > 0.0 {
        0.0
      } else {
        l
      }
    }
  }
}

///|
/// A test result whose null is continuous and whose two-sided p-value
/// doubles the smaller tail.
fn continuous_result(
  statistic : Double,
  null_dist : &@dist.Univariate,
  alternative : Alternative,
  df~ : Double?,
  df2~ : Double? = None,
  estimate~ : Double? = None,
  null_value~ : Double? = None,
  interval~ : Interval? = None,
  effect~ : Effect? = None,
  p_kind~ : PValueKind,
  procedure~ : String,
) -> TestResult {
  let nm = NullModel::new(observed=statistic, null_dist~, p_kind~, procedure~)
  {
    statistic,
    p_value: nm.p_value(alternative),
    log_p_value: log_tail(null_dist, statistic, alternative),
    df,
    df2,
    alternative,
    estimate,
    null_value,
    interval,
    effect,
    p_kind,
    null_model: Some(nm),
    procedure,
  }
}

///|
/// A test result on a discrete null: the p-value comes from the null
/// model, so the result and the model it keeps cannot disagree.
fn discrete_result(
  statistic : Double,
  nm : NullModel,
  alternative : Alternative,
  estimate~ : Double? = None,
  null_value~ : Double? = None,
  interval~ : Interval? = None,
  effect~ : Effect? = None,
) -> TestResult {
  let p = nm.p_value(alternative)
  {
    statistic,
    p_value: p,
    log_p_value: @math.ln(p),
    df: None,
    df2: None,
    alternative,
    estimate,
    null_value,
    interval,
    effect,
    p_kind: nm.p_kind,
    null_model: Some(nm),
    procedure: nm.procedure,
  }
}

///|
/// Combined ranks of several samples (average ties), the per-sample rank
/// sums, and the tie term Σ (t³ - t) over tied groups.
fn combined_ranks(
  groups : ArrayView[ArrayView[Double]],
) -> (FixedArray[Double], Double) {
  let all = []
  for g in groups {
    for v in g {
      all.push(v)
    }
  }
  let r = @describe.ranks(all[:], Average)
  let sums = FixedArray::make(groups.length(), 0.0)
  let mut i = 0
  for gi, g in groups {
    for _ in g {
      sums[gi] = sums[gi] + r[i]
      i = i + 1
    }
  }
  // tie term from the sorted values
  let sorted = FixedArray::from_iter(all.iter())
  sorted.sort()
  let mut ties = 0.0
  let mut a = 0
  while a < sorted.length() {
    let mut b = a
    while b + 1 < sorted.length() && sorted[b + 1] == sorted[a] {
      b = b + 1
    }
    let t = (b - a + 1).to_double()
    ties = ties + t * t * t - t
    a = b + 1
  }
  (sums, ties)
}
