///|
/// The tests, by family. Each takes the data shape from `data`, every
/// convention from `conventions.mbt`, and the assumptions it relies on as
/// certificates. Results are checked against SciPy by
/// infer_oracle_test.mbt.

// ---------------------------------------------------------------------------
// Location: t tests.
// ---------------------------------------------------------------------------

///|
fn t_from(
  est : Double,
  se : Double,
  df : Double,
  null_value : Double,
  alternative : Alternative,
  level : Double,
  effect : Effect?,
  procedure : String,
) -> TestResult raise QuestionError {
  if !(se > 0.0) {
    raise Degenerate(reason="zero standard error")
  }
  let t = @dist.StudentT::new(df) catch {
    _ => raise Degenerate(reason="degrees of freedom must be positive")
  }
  let null_dist = t as &@dist.Univariate
  continuous_result(
    (est - null_value) / se,
    null_dist,
    alternative,
    df=Some(df),
    estimate=Some(est),
    null_value=Some(null_value),
    interval=Some(
      pivot_interval(null_dist, est, se, level, alternative, procedure),
    ),
    effect~,
    p_kind=ExactNull,
    procedure~,
  )
}

///|
pub fn t_test_one(
  x : ArrayView[Double],
  mu~ : Double,
  alternative~ : Alternative,
  level~ : Double,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore(normality)
  check_level(level)
  let (n, m, v) = summarize(x)
  if n < 2 {
    raise TooFewPoints(needed=2, got=n)
  }
  let sd = v.sqrt()
  let effect : Effect = { kind: CohenD, value: (m - mu) / sd, interval: None, }
  t_from(
    m,
    sd / n.to_double().sqrt(),
    (n - 1).to_double(),
    mu,
    alternative,
    level,
    Some(effect),
    "one-sample t test",
  )
}

///|
fn two_sample_t(
  n1 : Int64,
  m1 : Double,
  v1 : Double,
  n2 : Int64,
  m2 : Double,
  v2 : Double,
  variances : Variances,
  alternative : Alternative,
  level : Double,
) -> TestResult raise QuestionError {
  check_level(level)
  if n1 < 2 || n2 < 2 {
    raise TooFewPoints(needed=2, got=(if n1 < n2 { n1 } else { n2 }).to_int())
  }
  let pooled = ((n1 - 1).to_double() * v1 + (n2 - 1).to_double() * v2) /
    (n1 + n2 - 2).to_double()
  let effect : Effect = {
    kind: CohenD,
    value: (m1 - m2) / pooled.sqrt(),
    interval: None,
  }
  match variances {
    Pooled(_) => {
      let se = (pooled * (1.0 / n1.to_double() + 1.0 / n2.to_double())).sqrt()
      t_from(
        m1 - m2,
        se,
        (n1 + n2 - 2).to_double(),
        0.0,
        alternative,
        level,
        Some(effect),
        "two-sample t test (pooled variance)",
      )
    }
    Welch => {
      let a = v1 / n1.to_double()
      let b = v2 / n2.to_double()
      let df = (a + b) *
        (a + b) /
        (a * a / (n1 - 1).to_double() + b * b / (n2 - 1).to_double())
      t_from(
        m1 - m2,
        (a + b).sqrt(),
        df,
        0.0,
        alternative,
        level,
        Some(effect),
        "Welch two-sample t test",
      )
    }
  }
}

///|
/// Two independent samples. `Pooled` carries the equal-variance
/// certificate, so the pooled test cannot be asked for without it.
pub fn t_test_two(
  x : ArrayView[Double],
  y : ArrayView[Double],
  variances~ : Variances,
  alternative~ : Alternative,
  level~ : Double,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore(normality)
  let (n1, m1, v1) = summarize(x)
  let (n2, m2, v2) = summarize(y)
  two_sample_t(
    n1.to_int64(),
    m1,
    v1,
    n2.to_int64(),
    m2,
    v2,
    variances,
    alternative,
    level,
  )
}

///|
/// The one-sample t test on the differences x - y of the pairs.
pub fn t_test_paired(
  xy : @data.Bivariate,
  alternative~ : Alternative,
  level~ : Double,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  let d = differences(xy)
  let r = t_test_one(d[:], mu=0.0, alternative~, level~, normality~)
  { ..r, procedure: "paired t test", }
}

///|
/// From published summaries only.
pub fn t_test_summary(
  x : @data.Summary,
  y : @data.Summary,
  variances~ : Variances,
  alternative~ : Alternative,
  level~ : Double,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore(normality)
  two_sample_t(
    x.n,
    x.mean,
    x.sd * x.sd,
    y.n,
    y.mean,
    y.sd * y.sd,
    variances,
    alternative,
    level,
  )
}

// ---------------------------------------------------------------------------
// Proportions and counts.
// ---------------------------------------------------------------------------

///|
/// Clopper–Pearson interval for k successes in n trials.
fn clopper_pearson(
  k : Int64,
  n : Int64,
  level : Double,
  alternative : Alternative,
) -> Interval raise QuestionError {
  let alpha = 1.0 - level
  let kd = k.to_double()
  let nd = n.to_double()
  let lower = (tail : Double) => {
    if k == 0L {
      0.0
    } else {
      let b = @dist.Beta::new(a=kd, b=nd - kd + 1.0) catch {
        _ => raise Degenerate(reason="bad counts")
      }
      @dist.Quantile::quantile(b, tail)
    }
  }
  let upper = (tail : Double) => {
    if k == n {
      1.0
    } else {
      let b = @dist.Beta::new(a=kd + 1.0, b=nd - kd) catch {
        _ => raise Degenerate(reason="bad counts")
      }
      @dist.Quantile::isf(b, tail)
    }
  }
  let procedure = "Clopper-Pearson"
  match alternative {
    TwoSided =>
      {
        lo: lower(0.5 * alpha),
        hi: upper(0.5 * alpha),
        level,
        kind: Confidence,
        procedure,
      }
    Greater =>
      { lo: lower(alpha), hi: 1.0, level, kind: Confidence, procedure, }
    Less => { lo: 0.0, hi: upper(alpha), level, kind: Confidence, procedure, }
  }
}

///|
/// Exact binomial test.
pub fn binomial_test(
  x : @data.Proportion,
  p0~ : Double,
  alternative~ : Alternative,
  two_sided~ : TwoSidedRule,
  level~ : Double,
) -> TestResult raise QuestionError {
  check_level(level)
  let n = x.trials
  let k = x.successes
  if n <= 0L || k < 0L || k > n {
    raise InvalidArgument(
      name="x",
      reason="need 0 <= successes <= trials, trials > 0",
    )
  }
  let d = @dist.Binomial::new(n~, p=p0) catch {
    _ => raise InvalidArgument(name="p0", reason="must be in [0, 1]")
  }
  let nm = NullModel::discrete(
    observed=k,
    null_dist=d,
    two_sided~,
    p_kind=ExactNull,
    procedure="exact binomial test",
  )
  discrete_result(
    k.to_double(),
    nm,
    alternative,
    estimate=Some(k.to_double() / n.to_double()),
    null_value=Some(p0),
    interval=Some(clopper_pearson(k, n, level, alternative)),
  )
}

///|
/// Score (z) test for one or several proportions.
pub fn proportion_test(
  xs : ArrayView[@data.Proportion],
  p0? : Double,
  continuity~ : Bool,
  alternative~ : Alternative,
  level~ : Double,
) -> TestResult raise QuestionError {
  ignore((xs, p0, continuity, alternative, level))
  raise NotImplemented(what="proportion_test")
}

///|
pub fn chi_square_goodness_of_fit(
  counts : @data.CategoryCounts,
  fitted_params~ : Int32,
) -> TestResult raise QuestionError {
  let k = counts.counts.length()
  if k < 2 {
    raise TooFewPoints(needed=2, got=k)
  }
  let mut total = 0.0
  for c in counts.counts {
    total = total + c.to_double()
  }
  let mut stat = @num.CompensatedSum::zero()
  for i in 0..<k {
    let p = match counts.null_probs {
      Some(ps) => ps[i]
      None => 1.0 / k.to_double()
    }
    let e = total * p
    if !(e > 0.0) {
      raise Degenerate(reason="an expected count is zero")
    }
    let o = counts.counts[i].to_double()
    stat = stat.add((o - e) * (o - e) / e)
  }
  let df = (k - 1 - fitted_params).to_double()
  let chi = @dist.ChiSquared::new(df) catch {
    _ => raise Degenerate(reason="no degrees of freedom left")
  }
  continuous_result(
    stat.value(),
    chi,
    Greater,
    df=Some(df),
    p_kind=AsymptoticNull,
    procedure="chi-square goodness-of-fit test",
  )
}

///|
/// Pearson's test of independence. With `continuity` and one degree of
/// freedom, each observed count moves towards its expectation by up to 1/2
/// (Yates), as R and SciPy do.
pub fn chi_square_independence(
  table : @data.Table,
  continuity~ : Bool,
) -> TestResult raise QuestionError {
  let r = table.rows
  let c = table.cols
  if r < 2 || c < 2 {
    raise TooFewPoints(needed=2, got=if r < c { r } else { c })
  }
  let total = table.total().to_double()
  let rows = FixedArray::make(r, 0.0)
  let cols = FixedArray::make(c, 0.0)
  for i in 0..<r {
    for j in 0..<c {
      let v = table.at(i, j).to_double()
      rows[i] = rows[i] + v
      cols[j] = cols[j] + v
    }
  }
  let df = ((r - 1) * (c - 1)).to_double()
  let mut stat = @num.CompensatedSum::zero()
  for i in 0..<r {
    for j in 0..<c {
      let e = rows[i] * cols[j] / total
      if !(e > 0.0) {
        raise Degenerate(reason="a row or column total is zero")
      }
      let mut o = table.at(i, j).to_double()
      if continuity && df == 1.0 {
        let diff = e - o
        let m = if diff.abs() < 0.5 { diff.abs() } else { 0.5 }
        o = o + (if diff > 0.0 { m } else if diff < 0.0 { -m } else { 0.0 })
      }
      stat = stat.add((o - e) * (o - e) / e)
    }
  }
  let chi = @dist.ChiSquared::new(df) catch {
    _ => raise Degenerate(reason="no degrees of freedom")
  }
  continuous_result(
    stat.value(),
    chi,
    Greater,
    df=Some(df),
    p_kind=AsymptoticNull,
    procedure="chi-square test of independence",
  )
}

///|
/// Fisher's exact test. 2 × 2 tables: hypergeometric null on the top-left
/// cell, estimate the sample odds ratio. Larger tables need the network
/// algorithm, not yet written.
pub fn fisher_exact(
  table : @data.Table,
  alternative~ : Alternative,
  two_sided~ : TwoSidedRule,
) -> TestResult raise QuestionError {
  if table.rows != 2 || table.cols != 2 {
    raise NotImplemented(what="Fisher's exact test beyond 2 x 2")
  }
  let a = table.at(0, 0)
  let b = table.at(0, 1)
  let c = table.at(1, 0)
  let d = table.at(1, 1)
  let m = a + b + c + d
  let h = @dist.Hypergeometric::new(population=m, successes=a + b, draws=a + c) catch {
    _ => raise Degenerate(reason="bad table")
  }
  let nm = NullModel::discrete(
    observed=a,
    null_dist=h,
    two_sided~,
    p_kind=ExactNull,
    procedure="Fisher's exact test",
  )
  let odds = a.to_double() * d.to_double() / (b.to_double() * c.to_double())
  discrete_result(
    odds,
    nm,
    alternative,
    estimate=Some(odds),
    null_value=Some(1.0),
    effect=Some({ kind: OddsRatio, value: odds, interval: None, }),
  )
}

///|
pub fn mcnemar(
  table : @data.Table,
  exactness~ : Exactness,
) -> TestResult raise QuestionError {
  ignore((table, exactness))
  raise NotImplemented(what="mcnemar")
}

// ---------------------------------------------------------------------------
// Several groups.
// ---------------------------------------------------------------------------

///|
fn anova_core(
  groups : ArrayView[ArrayView[Double]],
  pooled~ : Bool,
  procedure : String,
) -> TestResult raise QuestionError {
  let k = groups.length()
  if k < 2 {
    raise TooFewPoints(needed=2, got=k)
  }
  let ns = FixedArray::make(k, 0)
  let ms = FixedArray::make(k, 0.0)
  let vs = FixedArray::make(k, 0.0)
  let mut big_n = 0
  for i, g in groups {
    let (n, m, v) = summarize(g)
    if n < 2 {
      raise TooFewPoints(needed=2, got=n)
    }
    ns[i] = n
    ms[i] = m
    vs[i] = v
    big_n = big_n + n
  }
  match pooled {
    true => {
      let mut grand = 0.0
      for i in 0..<k {
        grand = grand + ns[i].to_double() * ms[i]
      }
      grand = grand / big_n.to_double()
      let mut ssb = 0.0
      let mut ssw = 0.0
      for i in 0..<k {
        ssb = ssb + ns[i].to_double() * (ms[i] - grand) * (ms[i] - grand)
        ssw = ssw + (ns[i] - 1).to_double() * vs[i]
      }
      let df1 = (k - 1).to_double()
      let df2 = (big_n - k).to_double()
      if !(ssw > 0.0) {
        raise Degenerate(reason="no variance within groups")
      }
      let f = @dist.FisherF::new(dfn=df1, dfd=df2) catch {
        _ => raise Degenerate(reason="too few observations")
      }
      continuous_result(
        ssb / df1 / (ssw / df2),
        f,
        Greater,
        df=Some(df1),
        df2=Some(df2),
        effect=Some({
          kind: EtaSquared,
          value: ssb / (ssb + ssw),
          interval: None,
        }),
        p_kind=ExactNull,
        procedure~,
      )
    }
    false => {
      let ws = FixedArray::makei(k, i => ns[i].to_double() / vs[i])
      let mut wsum = 0.0
      let mut wm = 0.0
      for i in 0..<k {
        wsum = wsum + ws[i]
        wm = wm + ws[i] * ms[i]
      }
      wm = wm / wsum
      let mut a = 0.0
      let mut lam = 0.0
      for i in 0..<k {
        a = a + ws[i] * (ms[i] - wm) * (ms[i] - wm)
        let r = 1.0 - ws[i] / wsum
        lam = lam + r * r / (ns[i] - 1).to_double()
      }
      let kd = k.to_double()
      let num = a / (kd - 1.0)
      let den = 1.0 + 2.0 * (kd - 2.0) / (kd * kd - 1.0) * lam
      let df1 = kd - 1.0
      let df2 = (kd * kd - 1.0) / (3.0 * lam)
      let f = @dist.FisherF::new(dfn=df1, dfd=df2) catch {
        _ => raise Degenerate(reason="too few observations")
      }
      continuous_result(
        num / den,
        f,
        Greater,
        df=Some(df1),
        df2=Some(df2),
        p_kind=AsymptoticNull,
        procedure="Welch's ANOVA",
      )
    }
  }
}

///|
/// One-way ANOVA; `Welch` gives Welch's heteroscedastic version, and the
/// classical one comes with its equal-variance certificate in `Pooled`.
pub fn anova_one_way(
  groups : ArrayView[ArrayView[Double]],
  variances~ : Variances,
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore(normality)
  anova_core(groups, pooled=variances is Pooled(_), "one-way ANOVA")
}

///|
/// Levene's original uses the mean; Brown–Forsythe the median.
#derive(Eq, Debug)
pub(all) enum LeveneCenter {
  Mean
  Median
  Trimmed(Double)
}

///|
/// ANOVA on the absolute deviations from each group's centre.
pub fn levene(
  groups : ArrayView[ArrayView[Double]],
  center~ : LeveneCenter,
) -> TestResult raise QuestionError {
  let devs : Array[ArrayView[Double]] = []
  for g in groups {
    let c = match center {
      Mean => summarize(g).1
      Median =>
        @describe.quantile(g, 0.5, Type7) catch {
          _ => raise Degenerate(reason="empty group")
        }
      Trimmed(p) =>
        @describe.trimmed_mean(g, p) catch {
          _ => raise Degenerate(reason="empty group")
        }
    }
    devs.push(FixedArray::from_iter(g.iter().map(v => (v - c).abs()))[:])
  }
  let r = anova_core(devs[:], pooled=true, "Levene's test")
  { ..r, effect: None, }
}

///|
pub fn bartlett(
  groups : ArrayView[ArrayView[Double]],
  normality~ : Assumed[Normality],
) -> TestResult raise QuestionError {
  ignore(normality)
  let k = groups.length()
  if k < 2 {
    raise TooFewPoints(needed=2, got=k)
  }
  let mut big_n = 0
  let mut sp = 0.0
  let mut sum_log = 0.0
  let mut sum_inv = 0.0
  for g in groups {
    let (n, _, v) = summarize(g)
    if n < 2 || !(v > 0.0) {
      raise Degenerate(reason="every group needs two values and some spread")
    }
    big_n = big_n + n
    sp = sp + (n - 1).to_double() * v
    sum_log = sum_log + (n - 1).to_double() * @math.ln(v)
    sum_inv = sum_inv + 1.0 / (n - 1).to_double()
  }
  let nk = (big_n - k).to_double()
  sp = sp / nk
  let kd = k.to_double()
  let stat = (nk * @math.ln(sp) - sum_log) /
    (1.0 + (sum_inv - 1.0 / nk) / (3.0 * (kd - 1.0)))
  let chi = @dist.ChiSquared::new(kd - 1.0) catch {
    _ => raise Degenerate(reason="need two groups")
  }
  continuous_result(
    stat,
    chi,
    Greater,
    df=Some(kd - 1.0),
    p_kind=AsymptoticNull,
    procedure="Bartlett's test",
  )
}

///|
/// Pairwise comparisons after ANOVA.
#derive(Eq, Debug)
pub(all) enum PostHoc {
  TukeyHsd
  Dunnett(control~ : Int32)
  /// Pairwise t tests with a p-value adjustment.
  Pairwise(Adjust)
}

///|
/// One `Comparison` per pair of groups, each naming its pair.
pub fn post_hoc(
  groups : ArrayView[ArrayView[Double]],
  how : PostHoc,
  level~ : Double,
) -> Array[Comparison] raise QuestionError {
  ignore((groups, how, level))
  raise NotImplemented(
    what="post-hoc comparisons (they need the studentized range)",
  )
}

// ---------------------------------------------------------------------------
// Rank tests.
// ---------------------------------------------------------------------------

///|
/// The exact null distribution of the Mann–Whitney U for sample sizes m
/// and n without ties: the coefficients of the Gaussian binomial
/// [m+n choose m]_q, built by multiplying by (1 - q^(n+i)) and dividing by
/// (1 - q^i) for i = 1..m.
fn mann_whitney_null(
  m : Int32,
  n : Int32,
) -> @dist.FiniteDiscrete raise QuestionError {
  let deg = m * n
  // after multiplying by (1 - q^(n+i)) and before dividing by (1 - q^i) the
  // degree is i(n + 1), up to m more than the final mn
  let c = FixedArray::make(deg + m + 1, 0.0)
  c[0] = 1.0
  let mut top = 0
  for i in 1..<=m {
    // multiply by (1 - q^(n+i)); the degree grows by n + i
    let shift = n + i
    let new_top = top + shift
    for j = new_top; j >= shift; j = j - 1 {
      c[j] = c[j] - c[j - shift]
    }
    top = new_top
    // divide by (1 - q^i): g[j] = f[j] + g[j - i]
    for j in i..<=top {
      c[j] = c[j] + c[j - i]
    }
    top = top - i
  }
  let mut total = 0.0
  for j in 0..<=deg {
    total = total + c[j]
  }
  let probs = FixedArray::makei(deg + 1, j => c[j] / total)
  @dist.FiniteDiscrete::new(0, probs[:]) catch {
    _ => raise Degenerate(reason="exact null failed")
  }
}

///|
/// The Mann–Whitney U test; the statistic is U for x. Exact when asked and
/// possible (no ties); otherwise the normal approximation with the tie
/// correction and, if asked, the continuity correction, as SciPy computes it.
pub fn mann_whitney(
  x : ArrayView[Double],
  y : ArrayView[Double],
  alternative~ : Alternative,
  exactness~ : Exactness,
  continuity~ : Bool,
) -> TestResult raise QuestionError {
  let n1 = x.length()
  let n2 = y.length()
  if n1 < 1 || n2 < 1 {
    raise TooFewPoints(needed=1, got=0)
  }
  let (sums, ties) = combined_ranks([x, y][:])
  let u1 = sums[0] - (n1 * (n1 + 1)).to_double() / 2.0
  let nn = (n1 * n2).to_double()
  let effect : Effect = {
    kind: RankBiserial,
    value: 2.0 * u1 / nn - 1.0,
    interval: None,
  }
  let exact = match exactness {
    Exact => {
      if ties > 0.0 {
        raise NotImplemented(what="exact Mann-Whitney null with ties")
      }
      true
    }
    Asymptotic => false
    ExactUpTo(k) => ties == 0.0 && n1 <= k && n2 <= k
  }
  if exact {
    let nm = NullModel::new(
      observed=u1,
      null_dist=mann_whitney_null(n1, n2),
      p_kind=ExactNull,
      procedure="Mann-Whitney U test",
    )
    return discrete_result(u1, nm, alternative, effect=Some(effect))
  }
  let n = (n1 + n2).to_double()
  let mu = nn / 2.0
  let sigma = (nn / 12.0 * (n + 1.0 - ties / (n * (n - 1.0)))).sqrt()
  if !(sigma > 0.0) {
    raise Degenerate(reason="all values tied")
  }
  let cc = if continuity { 0.5 } else { 0.0 }
  // SciPy: the upper tail of the larger U for two-sided, of U for x for
  // greater, of U for y for less
  let u = match alternative {
    TwoSided => if u1 > nn - u1 { u1 } else { nn - u1 }
    Greater => u1
    Less => nn - u1
  }
  let z = (u - mu - cc) / sigma
  let normal = @dist.Normal::new() catch { _ => panic() }
  let p1 = @dist.Cdf::sf(normal, z)
  let p = match alternative {
    TwoSided => if 2.0 * p1 > 1.0 { 1.0 } else { 2.0 * p1 }
    _ => p1
  }
  {
    statistic: u1,
    p_value: p,
    log_p_value: @math.ln(p),
    df: None,
    df2: None,
    alternative,
    estimate: None,
    null_value: None,
    interval: None,
    effect: Some(effect),
    p_kind: AsymptoticNull,
    null_model: None,
    procedure: "Mann-Whitney U test",
  }
}

///|
pub fn wilcoxon_signed_rank(
  x : ArrayView[Double],
  y? : ArrayView[Double],
  alternative~ : Alternative,
  zeros~ : Zeros,
  exactness~ : Exactness,
  continuity~ : Bool,
) -> TestResult raise QuestionError {
  ignore((x, y, alternative, zeros, exactness, continuity))
  raise NotImplemented(what="wilcoxon_signed_rank")
}

///|
/// Kruskal–Wallis H with the tie correction, against χ² with k - 1 df.
pub fn kruskal_wallis(
  groups : ArrayView[ArrayView[Double]],
) -> TestResult raise QuestionError {
  let k = groups.length()
  if k < 2 {
    raise TooFewPoints(needed=2, got=k)
  }
  let (sums, ties) = combined_ranks(groups)
  let mut n = 0
  let mut acc = 0.0
  for i, g in groups {
    if g.length() == 0 {
      raise TooFewPoints(needed=1, got=0)
    }
    n = n + g.length()
    acc = acc + sums[i] * sums[i] / g.length().to_double()
  }
  let nd = n.to_double()
  let h = 12.0 / (nd * (nd + 1.0)) * acc - 3.0 * (nd + 1.0)
  let correction = 1.0 - ties / (nd * nd * nd - nd)
  if !(correction > 0.0) {
    raise Degenerate(reason="all values tied")
  }
  let chi = @dist.ChiSquared::new((k - 1).to_double()) catch { _ => panic() }
  continuous_result(
    h / correction,
    chi,
    Greater,
    df=Some((k - 1).to_double()),
    p_kind=AsymptoticNull,
    procedure="Kruskal-Wallis test",
  )
}

///|
pub fn friedman(blocks : @data.Blocks) -> TestResult raise QuestionError {
  ignore(blocks)
  raise NotImplemented(what="friedman")
}

// ---------------------------------------------------------------------------
// Goodness of fit and normality.
// ---------------------------------------------------------------------------

///|
/// One-sample Kolmogorov–Smirnov against any distribution with a CDF.
pub fn ks_one(
  x : ArrayView[Double],
  null_dist : &@dist.Univariate,
  alternative~ : Alternative,
  exactness~ : Exactness,
) -> TestResult raise QuestionError {
  ignore((x, null_dist, alternative, exactness))
  raise NotImplemented(what="ks_one")
}

///|
pub fn ks_two(
  x : ArrayView[Double],
  y : ArrayView[Double],
  alternative~ : Alternative,
  exactness~ : Exactness,
) -> TestResult raise QuestionError {
  ignore((x, y, alternative, exactness))
  raise NotImplemented(what="ks_two")
}

///|
pub fn shapiro_wilk(x : ArrayView[Double]) -> TestResult raise QuestionError {
  ignore(x)
  raise NotImplemented(what="shapiro_wilk")
}

///|
pub fn anderson_darling(
  x : ArrayView[Double],
  null_dist : &@dist.Univariate,
) -> TestResult raise QuestionError {
  ignore((x, null_dist))
  raise NotImplemented(what="anderson_darling")
}

// ---------------------------------------------------------------------------
// Correlation and likelihood.
// ---------------------------------------------------------------------------

///|
#derive(Eq, Debug)
pub(all) enum CorrelationKind {
  Pearson
  Spearman
  KendallTauB
}

///|
/// Test of zero correlation against t with n - 2 df; for Pearson, the
/// Fisher z interval.
pub fn correlation_test(
  xy : @data.Bivariate,
  kind : CorrelationKind,
  alternative~ : Alternative,
  level~ : Double,
) -> TestResult raise QuestionError {
  check_level(level)
  let n = xy.x.length()
  if n < 3 {
    raise TooFewPoints(needed=3, got=n)
  }
  let (xs, ys) = match kind {
    Pearson => (xy.x, xy.y)
    Spearman =>
      (@describe.ranks(xy.x[:], Average), @describe.ranks(xy.y[:], Average))
    KendallTauB => raise NotImplemented(what="Kendall's tau test")
  }
  let mut c = @describe.Comoment::empty()
  for i in 0..<n {
    c = c.push(xs[i], ys[i])
  }
  let r = match c.correlation() {
    Some(r) => r
    None => raise Degenerate(reason="a variable is constant")
  }
  let df = (n - 2).to_double()
  let t = r * (df / ((1.0 - r) * (1.0 + r))).sqrt()
  let null_dist = @dist.StudentT::new(df) catch { _ => panic() }
  let interval = match kind {
    Pearson if n > 3 => {
      let z = @math.atanh(r)
      let se = 1.0 / (n - 3).to_double().sqrt()
      let normal = @dist.Normal::new() catch { _ => panic() }
      let iv = pivot_interval(normal, z, se, level, alternative, "Fisher z")
      Some({ ..iv, lo: @math.tanh(iv.lo), hi: @math.tanh(iv.hi), })
    }
    _ => None
  }
  let procedure = match kind {
    Pearson => "Pearson correlation test"
    _ => "Spearman correlation test"
  }
  let res = continuous_result(
    t,
    null_dist,
    alternative,
    df=Some(df),
    estimate=Some(r),
    null_value=Some(0.0),
    interval~,
    effect=Some({ kind: R, value: r, interval: None, }),
    p_kind=if kind == Pearson { ExactNull } else { AsymptoticNull },
    procedure~,
  )
  res
}

///|
/// Likelihood ratio test of nested models: 2(ℓ₁ − ℓ₀) against χ² with the
/// difference in free parameters.
pub fn likelihood_ratio(
  loglik_null~ : Double,
  loglik_alt~ : Double,
  df~ : Int32,
) -> TestResult raise QuestionError {
  let chi = @dist.ChiSquared::new(df.to_double()) catch {
    _ => raise InvalidArgument(name="df", reason="must be positive")
  }
  continuous_result(
    2.0 * (loglik_alt - loglik_null),
    chi,
    Greater,
    df=Some(df.to_double()),
    p_kind=AsymptoticNull,
    procedure="likelihood ratio test",
  )
}
