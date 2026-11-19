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
  let m : @infer.NullModel = {
    observed: 1.0,
    null_dist,
    two_sided: Doubling,
    p_kind: Exact,
    procedure: "toy permutation null",
  }
  assert_eq(m.p_value(Greater), 0.75)
  assert_eq(m.p_value(Less), 0.5)
  assert_eq(m.p_value(TwoSided), 1.0)
}
