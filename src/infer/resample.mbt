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
  check_level(level)
  if reps < 2 {
    raise TooFewPoints(needed=2, got=reps)
  }
  if interval == BCa {
    raise Unsupported(reason="BCa intervals are not implemented yet")
  }
  let draw = (x : ArrayView[Double]) => {
    FixedArray::makei(x.length(), _ => x[@num.uniform_below(rng, x.length().to_uint64()).to_int()])
  }
  let resample = () => match sample {
    One(x) => One(draw(x)[:])
    Two(x, y) => Two(draw(x)[:], draw(y)[:])
    Paired(x, y) => {
      let idx = FixedArray::makei(x.length(), _ => @num.uniform_below(rng, x.length().to_uint64()).to_int())
      Paired(
        FixedArray::makei(x.length(), i => x[idx[i]])[:],
        FixedArray::makei(x.length(), i => y[idx[i]])[:],
      )
    }
    Groups(gs) => Groups(gs.map(g => draw(g)[:]))
    Blocks(b) => Blocks(b)
  }
  if sample is Blocks(_) {
    raise Unsupported(reason="bootstrapping blocks is not implemented yet")
  }
  let theta = statistic(sample)
  let values = FixedArray::makei(reps, _ => statistic(resample()))
  let dist = @dist.Empirical::new(values[:], quantile_method=Type7) catch {
    _ => raise Degenerate(reason="the statistic was NaN on a resample")
  }
  let alpha = 1.0 - level
  let qlo = @dist.Quantile::quantile(dist, 0.5 * alpha)
  let qhi = @dist.Quantile::quantile(dist, 1.0 - 0.5 * alpha)
  let (lo, hi, name) = match interval {
    Percentile => (qlo, qhi, "bootstrap percentile")
    Basic => (2.0 * theta - qhi, 2.0 * theta - qlo, "bootstrap basic")
    BCa => (0.0, 0.0, "")
  }
  let (_, _, v) = summarize(values[:])
  {
    estimate: theta,
    interval: { lo, hi, level, kind: Confidence, procedure: name },
    se: v.sqrt(),
    distribution: dist,
    reps,
  }
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
  let observed = statistic(sample)
  let values : Array[Double] = []
  match (scheme, sample) {
    (Independent, Two(x, y)) => {
      let pooled = FixedArray::makei(x.length() + y.length(), i => if i < x.length() {
        x[i]
      } else {
        y[i - x.length()]
      })
      let n1 = x.length()
      match reps {
        None => {
          // every split of the pooled values into sizes n1 and n2
          let n = pooled.length()
          let pick = FixedArray::make(n1, 0)
          let mut count = 0
          fn visit(start : Int, depth : Int) -> Unit raise QuestionError {
            if depth == n1 {
              count = count + 1
              if count > 1000000 {
                raise Unsupported(reason="more than 10^6 splits; pass reps")
              }
              let chosen = FixedArray::make(n, false)
              for i in 0..<n1 {
                chosen[pick[i]] = true
              }
              let a = []
              let b = []
              for i in 0..<n {
                if chosen[i] {
                  a.push(pooled[i])
                } else {
                  b.push(pooled[i])
                }
              }
              values.push(statistic(Two(a[:], b[:])))
              return
            }
            for i in start..<(n - (n1 - depth) + 1) {
              pick[depth] = i
              visit(i + 1, depth + 1)
            }
          }
          visit(0, 0)
        }
        Some(r) =>
          for _ in 0..<r {
            @num.shuffle(rng, pooled)
            values.push(statistic(Two(pooled[0:n1], pooled[n1:])))
          }
      }
    }
    (WithinPairs, Paired(x, y)) => {
      let n = x.length()
      let a = FixedArray::from_iter(x.iter())
      let b = FixedArray::from_iter(y.iter())
      match reps {
        None => {
          if n > 20 {
            raise Unsupported(reason="more than 2^20 swaps; pass reps")
          }
          for mask in 0..<(1 << n) {
            for i in 0..<n {
              let swap = ((mask >> i) & 1) == 1
              a[i] = if swap { y[i] } else { x[i] }
              b[i] = if swap { x[i] } else { y[i] }
            }
            values.push(statistic(Paired(a[:], b[:])))
          }
        }
        Some(r) =>
          for _ in 0..<r {
            for i in 0..<n {
              let swap = (@num.uniform_below(rng, 2UL) == 1UL)
              a[i] = if swap { y[i] } else { x[i] }
              b[i] = if swap { x[i] } else { y[i] }
            }
            values.push(statistic(Paired(a[:], b[:])))
          }
      }
    }
    (Pairings, Paired(x, y)) => {
      let b = FixedArray::from_iter(y.iter())
      match reps {
        None => raise Unsupported(reason="exact pairings need n! arrangements; pass reps")
        Some(r) =>
          for _ in 0..<r {
            @num.shuffle(rng, b)
            values.push(statistic(Paired(x, b[:])))
          }
      }
    }
    _ => raise ShapeMismatch(estimand=MeanDifference)
  }
  // values within the tolerance of the observed statistic count as equal
  let tol = ties.relative * observed.abs()
  for i in 0..<values.length() {
    if (values[i] - observed).abs() <= tol {
      values[i] = observed
    }
  }
  let p_kind = match reps {
    None => PValueKind::Exact
    Some(r) => {
      // including the observed statistic gives the (k + 1)/(N + 1) p-value
      values.push(observed)
      MonteCarlo(draws=r, plus_one=true)
    }
  }
  let null_dist = @dist.Empirical::new(values[:], quantile_method=Type7) catch {
    _ => raise Degenerate(reason="the statistic was NaN on a permutation")
  }
  { observed, null_dist, two_sided: Doubling, p_kind, procedure: "permutation test" }
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
  let t0 = statistic(observed)
  let values = Array::makei(reps, _ => statistic(null_sampler(rng)))
  values.push(t0)
  let null_dist = @dist.Empirical::new(values[:], quantile_method=Type7) catch {
    _ => raise Degenerate(reason="the statistic was NaN on a draw")
  }
  {
    observed: t0,
    null_dist,
    two_sided: Doubling,
    p_kind: MonteCarlo(draws=reps, plus_one=true),
    procedure: "Monte Carlo test",
  }
}

///|
/// Leave-one-out estimates of bias and standard error.
pub fn jackknife(sample : Sample, statistic : Statistic) -> Estimate raise QuestionError {
  match sample {
    One(x) => {
      let n = x.length()
      if n < 2 {
        raise TooFewPoints(needed=2, got=n)
      }
      let theta = statistic(sample)
      let loo = FixedArray::makei(n, i => {
        let rest = FixedArray::makei(n - 1, j => if j < i { x[j] } else { x[j + 1] })
        statistic(One(rest[:]))
      })
      let (_, m, v) = summarize(loo[:])
      let nd = n.to_double()
      {
        value: theta - (nd - 1.0) * (m - theta),
        se: Some(((nd - 1.0) * (nd - 1.0) / nd * v).sqrt()),
        interval: None,
        n: n.to_int64(),
        procedure: "jackknife (bias-corrected)",
      }
    }
    _ => raise Unsupported(reason="jackknife for one sample only, for now")
  }
}
