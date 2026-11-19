///|
/// Shared pieces of the tests: two-pass summaries, critical values,
/// result assembly, and the minimum-likelihood two-sided p-value for
/// unimodal discrete nulls.

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
    raise Degenerate(reason="level must be in (0, 1)")
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
  let nm : NullModel = {
    observed: statistic,
    null_dist,
    two_sided: Doubling,
    p_kind,
    procedure,
  }
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
/// Two-sided p-value of a unimodal discrete null by minimum likelihood: the
/// total probability of outcomes no more likely than the observed k, with
/// the relative slack 1 + 1e-7 that R and SciPy use so that outcomes equal
/// in exact arithmetic are not split by rounding.
fn[D : @dist.Mass + @dist.IntCdf] min_likelihood_p(
  d : D,
  k : Int64,
  lo : Int64,
  hi : Int64,
  mode : Double,
) -> Double {
  let bound = @dist.Mass::pmf(d, k) * (1.0 + 1.0e-7)
  let p = if k.to_double() < mode {
    // the other tail starts at the first j above the mode with pmf <= bound
    let mut j = @math.ceil(mode).to_int64()
    while j <= hi && @dist.Mass::pmf(d, j) > bound {
      j = j + 1
    }
    @dist.IntCdf::cdf_int(d, k) + @dist.IntCdf::sf_int(d, j - 1)
  } else if k.to_double() > mode {
    let mut j = @math.floor(mode).to_int64()
    while j >= lo && @dist.Mass::pmf(d, j) > bound {
      j = j - 1
    }
    @dist.IntCdf::cdf_int(d, j) + @dist.IntCdf::sf_int(d, k - 1)
  } else {
    1.0
  }
  if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Tail p-values of a discrete null at k: (P(X >= k), P(X <= k)).
fn[D : @dist.IntCdf] discrete_tails(d : D, k : Int64) -> (Double, Double) {
  (@dist.IntCdf::sf_int(d, k - 1), @dist.IntCdf::cdf_int(d, k))
}

///|
/// Choose the p-value of a discrete null for the alternative and rule.
fn[D : @dist.Mass + @dist.IntCdf] discrete_p(
  d : D,
  k : Int64,
  lo : Int64,
  hi : Int64,
  mode : Double,
  alternative : Alternative,
  rule : TwoSidedRule,
) -> Double {
  let (upper, lower) = discrete_tails(d, k)
  match alternative {
    Greater => upper
    Less => lower
    TwoSided =>
      match rule {
        Doubling => {
          let p = 2.0 * (if upper < lower { upper } else { lower })
          if p > 1.0 {
            1.0
          } else {
            p
          }
        }
        MinLikelihood => min_likelihood_p(d, k, lo, hi, mode)
      }
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
