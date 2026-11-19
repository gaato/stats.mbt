///|
fn near(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol * b.abs().max(1.0)
}

///|
test "constraints round-trip through the unconstrained line" {
  let cases : Array[(@dist.Constraint, Double)] = [
    (Real, -3.5),
    (Positive, 0.02),
    (UnitInterval, 0.9),
    (Interval(lo=-1.0, hi=4.0), 3.0),
  ]
  for case in cases {
    let (c, x) = case
    assert_true(c.contains(x))
    assert_true(near(c.from_unconstrained(c.to_unconstrained(x)), x, 1.0e-14))
  }
  assert_true(!@dist.Constraint::Positive.contains(0.0))
}

///|
test "a family rebuilds its members from parameters" {
  let fam = @dist.Normal::family()
  let d = (fam.make)([1.5, 2.0])
  assert_eq((fam.params)(d), [1.5, 2.0])
  // every continuous family has its shape, then loc and scale
  let gamma = @dist.Gamma::family()
  assert_eq(gamma.specs.map(s => s.name), ["a", "loc", "scale"])
  // gamma by rate and by scale are different parameterisations of the
  // same distributions
  let g = (gamma.make)([3.0, 0.0, 0.5])
  assert_eq((@dist.Gamma::family_rate().params)(g), [3.0, 2.0])
  // fixed hyperparameters live in the family value
  let b = @dist.Binomial::family(n=10)
  assert_eq(b.dim(), 1)
  assert_eq(b.specs[0].name, "p")
}

///|
test "normal density in closed form" {
  let d = @dist.Normal::new(loc=1.0, scale=2.0)
  assert_true(near(@dist.Density::pdf(d, 1.0), 0.19947114020071635, 1.0e-15))
  // the log density stays finite where the density underflows
  assert_true(@dist.Density::pdf(d, 1.0e3) == 0.0)
  assert_true(
    near(@dist.Density::logpdf(d, 1.0e3), -124751.73708571376, 1.0e-15),
  )
}

///|
test "normal-normal conjugate update" {
  let prior = @dist.Normal::new(loc=0.0, scale=1.0)
  let c = @dist.Normal::known_variance_prior(1.0)
  let post = (c.update)(prior, [2.0, 2.0])
  // precision 1 + 2 = 3, mean (0 + 4) / 3
  assert_true(near(@dist.Moments::mean(post).unwrap(), 4.0 / 3.0, 1.0e-15))
  assert_true(near(@dist.Moments::variance(post).unwrap(), 1.0 / 3.0, 1.0e-15))
}
