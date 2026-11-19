///|
/// The Double interface of a discrete distribution does not pass through
/// an integer type, so it is right beyond every integer range.
test "discrete quantiles and tails beyond the integer ranges" {
  // median ceil(ln(1/2)/ln(1 - p)): past Int32
  let g = @dist.Geometric::new(1.0e-12)
  assert_eq(g.quantile(0.5), 693147180560.0)
  assert_eq(g.quantile_int(0.5), 693147180560L)
  assert_eq(g.cdf(693147180560.0) >= 0.5, true)
  assert_eq(g.cdf(693147180559.0) < 0.5, true)
  // past Int64: the Double answer is exact, the integer one saturates
  let tiny = @dist.Geometric::new(1.0e-25)
  assert_true(
    (tiny.quantile(0.5) - 6.931471805599453e24).abs() <=
    1.0e-15 * 6.931471805599453e24,
  )
  assert_eq(tiny.quantile_int(0.5), 9223372036854775807L)
  assert_true((tiny.cdf(6.931471805599453e24) - 0.5).abs() < 1.0e-15)
  assert_true((tiny.sf(6.931471805599453e24) - 0.5).abs() < 1.0e-15)
  // a Poisson mean of 1e10: the distribution function at the mean
  let p = @dist.Poisson::new(1.0e10)
  assert_true((p.cdf(1.0e10) - 0.5000026596152026).abs() < 1.0e-12)
  assert_eq(p.quantile(0.5), 1.0e10)
}

///|
/// The upper quantile of a discrete distribution is found on the upper
/// tail, so a probability such as 1e-30 keeps its meaning.
test "discrete upper quantiles on the upper tail" {
  let g = @dist.Geometric::new(0.5)
  // P(X > k) = 2^-k: the smallest k with 2^-k <= q
  assert_eq(g.isf(0.5), 1.0)
  assert_eq(g.isf(0.25), 2.0)
  assert_eq(g.isf(0.3), 2.0)
  assert_eq(g.isf(1.0e-30), 100.0)
  assert_eq(g.isf(1.0), 1.0)
  assert_eq(g.isf(0.0), 1.0 / 0.0)
  let b = @dist.Binomial::new(n=10, p=0.5)
  // P(X > 9) = 2^-10, P(X > 8) = 11 · 2^-10
  assert_eq(b.isf(0.001), 9.0)
  assert_eq(b.isf(0.0009), 10.0)
  assert_eq(b.isf(0.02), 8.0)
  assert_eq(b.isf(1.0), 0.0)
  let l = @dist.DiscreteLaplace::new(1.0)
  for q in [0.999, 0.9, 0.6, 0.4, 0.1, 1.0e-5] {
    let k = l.isf(q)
    assert_true(l.sf(k) <= q, msg="sf(isf(\{q}))")
    assert_true(l.sf(k - 1.0) > q, msg="sf(isf(\{q}) - 1)")
  }
}

///|
test "supports, modes and moments of the discrete reference distributions" {
  let b = @dist.Binomial::new(n=10, p=0.3)
  assert_eq((b.support_lower(), b.support_upper()), (0.0, 10.0))
  assert_eq(b.mode(), 3.0)
  // (n + 1) p an integer: two equal masses, the smaller point
  let tie = @dist.Binomial::new(n=9, p=0.5)
  assert_eq(tie.mode(), 4.0)
  assert_eq(tie.pmf(4), tie.pmf(5))
  let p = @dist.Poisson::new(3.0)
  assert_eq(p.mode(), 2.0)
  assert_eq(@dist.Poisson::new(3.5).mode(), 3.0)
  assert_eq(p.mean(), Some(3.0))
  assert_eq(p.support_upper(), 1.0 / 0.0)
  let h = @dist.Hypergeometric::new(population=50, successes=12, draws=20)
  assert_eq((h.support_lower(), h.support_upper()), (0.0, 12.0))
  assert_true((h.mean().unwrap() - 4.8).abs() < 1.0e-15)
  assert_true((h.variance().unwrap() - 2.233469387755102).abs() < 1.0e-14)
  let m = h.mode()
  assert_true(h.pmf(m.to_int64()) >= h.pmf(m.to_int64() - 1))
  assert_true(h.pmf(m.to_int64()) >= h.pmf(m.to_int64() + 1))
  // the upper bound of a discrete uniform is not included
  let die = @dist.DiscreteUniform::new(low=1, high=7)
  assert_eq((die.support_lower(), die.support_upper()), (1.0, 6.0))
  assert_eq(die.pmf(6), 1.0 / 6.0)
  assert_eq(die.pmf(7), 0.0)
}

///|
test "a table of probabilities on the integers" {
  let d = @dist.FiniteDiscrete::new(-1, [0.0, 0.25, 0.25, 0.5, 0.0])
  assert_eq(d.offset(), -1L)
  assert_eq(d.probs().to_owned(), [0.0, 0.25, 0.25, 0.5, 0.0])
  assert_eq((d.support_lower(), d.support_upper()), (0.0, 2.0))
  assert_eq(d.pmf(2), 0.5)
  assert_eq(d.cdf(1.5), 0.5)
  assert_eq(d.sf(0.0), 0.75)
  assert_eq(d.sf_at_least(2.0), 0.5)
  assert_eq(d.quantile(0.5), 1.0)
  assert_eq(d.quantile_int(0.6), 2L)
  // the smallest k with P(X > k) <= q
  assert_eq(d.isf(0.5), 1.0)
  assert_eq(d.isf(0.4), 2.0)
  assert_eq(d.isf(0.75), 0.0)
  // a point far outside the table needs no integer conversion
  assert_eq(d.cdf(1.0e300), 1.0)
  assert_eq(d.sf(-1.0e300), 1.0)
}

///|
/// Atoms off the integers: the exact null of a rank statistic with ties
/// takes half-integer values.
test "atoms anywhere on the line" {
  let a = @dist.Atoms::new([2.5, 0.5, 1.0, 2.5], [0.125, 0.25, 0.5, 0.125])
  assert_eq(a.values().to_owned(), [0.5, 1.0, 2.5])
  assert_eq(a.probs().to_owned(), [0.25, 0.5, 0.25])
  assert_eq(a.mass(1.0), 0.5)
  assert_eq(a.mass(1.5), 0.0)
  assert_eq(a.cdf(0.75), 0.25)
  assert_eq(a.cdf(1.0), 0.75)
  assert_eq(a.sf(1.0), 0.25)
  assert_eq(a.sf_at_least(1.0), 0.75)
  assert_eq(a.sf_at_least(3.0), 0.0)
  assert_eq(a.quantile(0.25), 0.5)
  assert_eq(a.quantile(0.26), 1.0)
  assert_eq(a.isf(0.25), 1.0)
  assert_eq(a.isf(0.2), 2.5)
  assert_eq(a.mean(), Some(1.25))
  assert_eq((a.support_lower(), a.support_upper()), (0.5, 2.5))
  let bad = try @dist.Atoms::new([1.0, 2.0], [0.5, 0.25]) catch {
    InvalidParam(reason~, ..) => reason == "must sum to 1"
    @dist.ParamError::.. => false
  } noraise {
    _ => false
  }
  assert_true(bad)
}

///|
test "atoms from an empirical distribution and from weights" {
  let e = @dist.Empirical::new([3.0, 1.0, 3.0, 2.0], quantile_method=Type7)
  assert_eq(e.values().to_owned(), [1.0, 2.0, 3.0, 3.0])
  assert_eq(e.length(), 4)
  assert_eq(e.mass(3.0), 0.5)
  assert_eq(e.mass(2.5), 0.0)
  assert_eq((e.support_lower(), e.support_upper()), (1.0, 3.0))
  let a = e.to_atoms()
  assert_eq(a.values().to_owned(), [1.0, 2.0, 3.0])
  assert_eq(a.probs().to_owned(), [0.25, 0.25, 0.5])
  assert_eq(a.cdf(2.0), e.cdf(2.0))
  assert_eq(a.sf_at_least(3.0), e.sf_at_least(3.0))
  // frequency weights 1, 3: the same as repeating the second value
  let w = @data.Weighted::new(
    [10.0, 20.0],
    @data.Weights::new([1.0, 3.0], kind=Frequency),
  )
  let f = @dist.Atoms::from_weighted(w)
  assert_eq(f.probs().to_owned(), [0.25, 0.75])
  assert_eq(f.mean(), Some(17.5))
}
