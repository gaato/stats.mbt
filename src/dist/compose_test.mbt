///|
fn close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol
}

///|
/// Zero-inflation as a mixture of two atoms: the distribution function
/// jumps by each weight, and the upper tail comes from the components.
test "mixture of point masses" {
  let m = @dist.Mixture::new([0.3, 0.7], [
    @dist.PointMass::new(0.0) as &@dist.Univariate,
    @dist.PointMass::new(1.0),
  ])
  assert_eq(@dist.Cdf::cdf(m, -0.5), 0.0)
  assert_eq(@dist.Cdf::cdf(m, 0.5), 0.3)
  assert_eq(@dist.Cdf::sf(m, 0.5), 0.7)
  assert_eq(@dist.Cdf::cdf(m, 1.0), 1.0)
}

///|
test "mixture weights are checked, not renormalised" {
  let bad = try {
    @dist.Mixture::new([0.3, 0.3], [
      @dist.PointMass::new(0.0) as &@dist.Univariate,
      @dist.PointMass::new(1.0),
    ])
    |> ignore
    false
  } catch {
    @dist.InvalidParam(name~, ..) => name == "weights"
    _ => false
  }
  assert_true(bad)
}

///|
test "affine and truncated compose" {
  let m = @dist.Mixture::new([0.25, 0.25, 0.5], [
    @dist.PointMass::new(0.0) as &@dist.Univariate,
    @dist.PointMass::new(1.0),
    @dist.PointMass::new(2.0),
  ])
  // 10 + 2X moves the atoms to 10, 12, 14.
  let a = @dist.Affine::new(m, loc=10.0, scale=2.0)
  assert_eq(@dist.Cdf::cdf(a, 11.0), 0.25)
  assert_eq(@dist.Cdf::cdf(a, 12.0), 0.5)
  // Conditioning on 0 < X <= 2 removes the atom at 0.
  let t = @dist.Truncated::new(m, lo=0.0, hi=2.0)
  assert_true(close(@dist.Cdf::cdf(t, 1.5), 1.0 / 3.0, 1.0e-15))
}

///|
/// The polar-form power must agree with the defining sum.
test "binomial characteristic function" {
  let n = 3
  let p = 0.3
  let b = @dist.Binomial::new(n, p)
  for t in [0.0, 0.7, 2.5, -4.0] {
    let (re, im) = @dist.Characteristic::cf(b, t)
    let coef = [1.0, 3.0, 3.0, 1.0]
    let mut sre = 0.0
    let mut sim = 0.0
    for k in 0..<4 {
      let w = coef[k] *
        @math.pow(p, k.to_double()) *
        @math.pow(1.0 - p, (n - k).to_double())
      sre = sre + w * @math.cos(t * k.to_double())
      sim = sim + w * @math.sin(t * k.to_double())
    }
    assert_true(close(re, sre, 1.0e-14))
    assert_true(close(im, sim, 1.0e-14))
  }
}

///|
test "normal characteristic and moment generating functions" {
  let d = @dist.Normal::new(mu=1.5, sigma=2.0)
  let (re, im) = @dist.Characteristic::cf(d, 0.4)
  // |φ(t)| = exp(-σ²t²/2), arg φ(t) = μt.
  assert_true(close(@math.hypot(re, im), @math.exp(-0.5 * 4.0 * 0.16), 1.0e-15))
  assert_true(close(@math.atan2(im, re), 1.5 * 0.4, 1.0e-15))
  assert_eq(@dist.Mgf::mgf(d, 0.0), Some(1.0))
}

///|
test "gamma moment generating function diverges at the rate" {
  let g = @dist.Gamma::with_rate(2.0, 3.0)
  assert_eq(@dist.Mgf::mgf(g, 3.0), None)
  assert_eq(@dist.Mgf::cgf(g, 0.0), Some(0.0))
  assert_true(@dist.Mgf::mgf(g, 1.0) is Some(_))
}

///|
test "cantor moments exist without a density" {
  let c = @dist.Cantor::new()
  assert_eq(@dist.Moments::mean(c), Some(0.5))
  assert_eq(@dist.Moments::variance(c), Some(0.125))
}
