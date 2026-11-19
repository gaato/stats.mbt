///|
test "questions are checked against the data shape" {
  let a = [1.0, 2.0, 3.0]
  let wrong = try @infer.Question::new(MeanDifference, One(a[:])) catch {
    ShapeMismatch(..) => true
    _ => false
  } noraise {
    _ => false
  }
  assert_true(wrong)
  let q = @infer.Question::new(MeanDifference, Two(a[:], a[:]))
  assert_eq(q.estimand, MeanDifference)
}

///|
/// The same p-value code serves a randomisation null: here a discrete
/// null with atoms at 0, 1 and 2.
test "p-values from a discrete null count the observed value" {
  let null_dist = @dist.Mixture::new([0.25, 0.25, 0.5], [
    @dist.PointMass::new(0.0) as &@dist.Univariate,
    @dist.PointMass::new(1.0),
    @dist.PointMass::new(2.0),
  ])
  let m = @infer.NullModel::new(
    observed=1.0,
    null_dist~,
    p_kind=ExactNull,
    procedure="toy permutation null",
  )
  assert_eq(m.p_value(Greater), 0.75)
  assert_eq(m.p_value(Less), 0.5)
  assert_eq(m.p_value(TwoSided), 1.0)
}

///|
/// A question whose answer is not one number is refused when asked.
test "questions about several groups are not questions about one value" {
  let a = [1.0, 2.0, 3.0]
  let b = [2.0, 3.0, 5.0]
  let groups = try @infer.Question::new(Mean, Groups([a[:], b[:]])) catch {
    ShapeMismatch(estimand~) => estimand == Mean
    _ => false
  } noraise {
    _ => false
  }
  assert_true(groups)
  let p = try @infer.Question::new(Quantile(1.5), One(a[:])) catch {
    InvalidArgument(name~, ..) => name == "p"
    _ => false
  } noraise {
    _ => false
  }
  assert_true(p)
}

///|
/// Pairs carry their equal length, so resampling them cannot index past
/// the shorter side.
test "paired samples are pairs" {
  let rng = @num.Xoshiro256::seed(1)
  let xy = @data.Bivariate::new([1.0, 2.0, 3.0, 4.0, 5.0], [
    2.0, 2.5, 3.5, 5.0, 5.5,
  ])
  let q = @infer.Question::new(MeanDifference, Paired(xy))
  let e = @infer.estimate(q, Bootstrap(reps=200), level=0.9, rng~)
  // the mean of the differences x - y: -3.5 / 5
  assert_true((e.value - -0.7).abs() < 1.0e-15)
  assert_eq(e.n, 5L)
}

///|
/// The same two-sided p-value from the test and from the null model it
/// keeps, under either rule.
test "a discrete null model answers every alternative" {
  let r = @infer.binomial_test(
    { successes: 3, trials: 20, },
    p0=0.4,
    alternative=TwoSided,
    two_sided=MinLikelihood,
    level=0.95,
  )
  let nm = r.null_model.unwrap()
  assert_eq(nm.p_value(TwoSided), r.p_value)
  assert_true(nm.p_value(Less) < nm.p_value(TwoSided))
  assert_eq(nm.two_sided, MinLikelihood)
  // 2^32 trials: counts are Int64 all the way to the distribution
  let big = @infer.binomial_test(
    { successes: 2147483648L, trials: 4294967296L, },
    p0=0.5,
    alternative=TwoSided,
    two_sided=Doubling,
    level=0.95,
  )
  assert_eq(big.p_value, 1.0)
  assert_eq(big.estimate, Some(0.5))
}

///|
/// Atoms off the integers, with the minimum-likelihood rule: the atoms no
/// more likely than the observed one.
test "a null model on atoms" {
  let null_dist = @dist.Atoms::new([0.5, 1.0, 1.5, 2.0], [0.1, 0.4, 0.3, 0.2])
  let m = @infer.NullModel::atoms(
    observed=2.0,
    null_dist~,
    two_sided=MinLikelihood,
    p_kind=ExactNull,
    procedure="toy rank null",
  )
  // atoms with mass <= 0.2: 0.5 and 2.0
  assert_true((m.p_value(TwoSided) - 0.3).abs() < 1.0e-15)
  assert_eq(m.p_value(Greater), 0.2)
  debug_inspect(
    m,
    content=(
      #|{
      #|  observed: 2,
      #|  two_sided: MinLikelihood,
      #|  p_kind: ExactNull,
      #|  procedure: "toy rank null",
      #|}
    ),
  )
}

///|
test "arguments outside their range are named" {
  let rng = @num.Xoshiro256::seed(1)
  let a = [1.0, 2.0, 3.0]
  let level = try
    @infer.bootstrap(
      One(a[:]),
      _ => 0.0,
      reps=10,
      level=1.5,
      interval=Percentile,
      rng~,
    )
  catch {
    InvalidArgument(name~, ..) => name == "level"
    _ => false
  } noraise {
    _ => false
  }
  assert_true(level)
  let reps = try
    @infer.bootstrap(
      One(a[:]),
      _ => 0.0,
      reps=1,
      level=0.9,
      interval=Percentile,
      rng~,
    )
  catch {
    InvalidArgument(name~, ..) => name == "reps"
    _ => false
  } noraise {
    _ => false
  }
  assert_true(reps)
  let scheme = try
    @infer.permutation_test(
      One(a[:]),
      _ => 0.0,
      scheme=Pairings,
      reps=Random(10),
      ties={ relative: 0.0, },
      rng~,
    )
  catch {
    SchemeMismatch(scheme~) => scheme == Pairings
    _ => false
  } noraise {
    _ => false
  }
  assert_true(scheme)
}

///|
test "results print" {
  let normality : @infer.Assumed[@infer.Normality] = @infer.Assumed::postulate()
  let r = @infer.t_test_one(
    [1.0, 2.0, 4.0],
    mu=0.0,
    alternative=TwoSided,
    level=0.95,
    normality~,
  )
  let text = @debug.to_string(r)
  assert_true(text.contains("one-sample t test"))
}
