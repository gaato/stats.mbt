///|
fn rel(got : Double, want : Double, tol : Double) -> Bool {
  (got - want).abs() <= tol * want.abs()
}

///|
/// Methods are reachable with dot syntax, parameters can be read back, and
/// a distribution can be printed.
test "a distribution is used like an ordinary value" {
  let d = @dist.Normal::new(loc=1.0, scale=2.0)
  assert_eq(d.cdf(1.0), 0.5)
  assert_eq(d.sf(1.0), 0.5)
  assert_eq(d.quantile(0.5), 1.0)
  assert_eq(d.mean(), Some(1.0))
  assert_eq(d.variance(), Some(4.0))
  assert_eq(d.mode(), 1.0)
  assert_eq(d.loc, 1.0)
  assert_eq(d.scale, 2.0)
  debug_inspect(d, content="{ loc: 1, scale: 2 }")
  let b = @dist.Beta::new(a=2.0, b=3.0)
  assert_eq((b.a, b.b), (2.0, 3.0))
  let t = @dist.StudentT::new(5.0)
  assert_eq(t.df, 5.0)
  assert_true(d == @dist.Normal::new(loc=1.0, scale=2.0))
}

///|
/// Every continuous family has a location and a scale: the exponential
/// with rate 2 is the exponential with scale 1/2.
test "location and scale on every continuous family" {
  let e = @dist.Exponential::new(scale=0.5)
  assert_true(rel(e.cdf(1.0), 0.8646647167633873, 1.0e-15))
  assert_true(rel(e.isf(1.0e-30), 34.538776394910684, 1.0e-15))
  assert_true(rel(e.logsf(1000.0), -2000.0, 1.0e-15))
  let w = @dist.Weibull::new(1.5, loc=1.0, scale=2.0)
  assert_true(rel(w.sf(4.0), 0.15927590849002138, 1.0e-14))
  assert_eq(w.cdf(0.5), 0.0)
  assert_eq(w.support_lower(), 1.0)
  let u = @dist.Uniform::new(loc=1.0, scale=3.0)
  assert_eq((u.support_lower(), u.support_upper()), (1.0, 4.0))
  assert_eq(u.pdf(2.0), 1.0 / 3.0)
  assert_eq(u.quantile(0.5), 2.5)
  // the scale must be positive and finite
  let bad = try @dist.Exponential::new(scale=0.0) catch {
    InvalidParam(name~, ..) => name == "scale"
  } noraise {
    _ => false
  }
  assert_true(bad)
}

///|
/// The upper tail of a shifted and scaled family is the family's own
/// upper tail: a t quantile at 1e-20 is not rounded through 1 - q, and
/// the normal log tails stay finite where the tails underflow.
test "tails survive location and scale" {
  let t = @dist.StudentT::new(5.0, loc=3.0, scale=2.0)
  assert_true(rel(t.isf(1.0e-20), 3.0 + 2.0 * 15683.925454365775, 1.0e-14))
  let z = @dist.Normal::new()
  assert_eq(z.sf(40.0), 0.0)
  assert_true(rel(z.logsf(40.0), -804.6084420137538, 1.0e-15))
  assert_true(rel(z.logcdf(-40.0), -804.6084420137538, 1.0e-15))
  assert_true(rel(z.logpdf(40.0), -800.9189385332047, 1.0e-15))
  let n = @dist.Normal::new(loc=1.0, scale=2.0)
  assert_true(rel(n.logsf(81.0), -804.6084420137538, 1.0e-15))
  assert_true(rel(n.logpdf(81.0), -800.9189385332047 - @math.ln(2.0), 1.0e-15))
  // a wrapper forwards every tail too
  let a = @dist.Affine::new(@dist.StudentT::new(5.0), loc=3.0, scale=2.0)
  assert_eq(a.isf(1.0e-20), t.isf(1.0e-20))
  let an = @dist.Affine::new(z, loc=1.0, scale=2.0)
  assert_eq(an.logsf(81.0), n.logsf(81.0))
  assert_eq(an.logpdf(81.0), n.logpdf(81.0))
  assert_eq(an.mean(), Some(1.0))
  assert_eq(an.variance(), Some(4.0))
}

///|
/// A wrapper over a distribution is a distribution: it goes wherever a
/// `&Univariate` is taken.
test "wrappers are univariate distributions" {
  let a = @dist.Affine::new(@dist.StudentT::new(5.0), loc=3.0, scale=2.0)
  let t = @dist.Truncated::new(@dist.Normal::new(), lo=0.0)
  let c = @dist.Censored::new(@dist.Normal::new(), lo=-1.0, hi=1.0)
  let parts : Array[&@dist.Univariate] = [a, t, c]
  let m = @dist.Mixture::new([0.2, 0.3, 0.5], parts)
  assert_true(m.cdf(0.0) > 0.0 && m.cdf(0.0) < 1.0)
}

///|
/// Truncation far in the upper tail: the mass 1.1e-19 is a difference of
/// upper tails, where a difference of distribution functions is 0.
test "truncation keeps the tail it cuts into" {
  let t = @dist.Truncated::new(@dist.Normal::new(), lo=9.0)
  assert_true(rel(t.cdf(9.5), 0.9907012024755889, 1.0e-13))
  assert_true(rel(t.sf(9.5), 0.00929879752441107, 1.0e-13))
  assert_true(rel(t.sf(12.0), 1.5740743947978652e-14, 1.0e-13))
  assert_true(rel(t.quantile(0.5), 9.075787065491832, 1.0e-13))
  assert_true(rel(t.isf(1.0e-10), 11.252265916368975, 1.0e-13))
  assert_eq(t.cdf(9.0), 0.0)
  assert_eq(t.pdf(8.0), 0.0)
  assert_eq((t.support_lower(), t.support_upper()), (9.0, 1.0 / 0.0))
  // both bounds, around the centre
  let m = @dist.Truncated::new(@dist.Normal::new(), lo=-1.0, hi=2.0)
  assert_true(rel(m.cdf(0.5), 0.6508804213366272, 1.0e-14))
  assert_true(rel(m.sf(0.5), 0.34911957866337284, 1.0e-14))
  assert_true(rel(m.quantile(m.cdf(0.5)), 0.5, 1.0e-13))
  assert_true(rel(m.isf(m.sf(0.5)), 0.5, 1.0e-13))
  // an interval the base gives no mass is refused
  let empty = try
    @dist.Truncated::new(@dist.Exponential::new(), lo=-3.0, hi=-1.0)
  catch {
    InvalidParam(reason~, ..) => reason == "interval has no mass"
  } noraise {
    _ => false
  }
  assert_true(empty)
}

///|
test "censoring moves the tails onto the limits" {
  let c = @dist.Censored::new(@dist.Normal::new(), lo=-1.0, hi=1.0)
  let z = @dist.Normal::new()
  assert_eq(c.cdf(-1.0), z.cdf(-1.0))
  assert_eq(c.cdf(-1.5), 0.0)
  assert_eq(c.sf(1.0), 0.0)
  assert_eq(c.sf_at_least(1.0), z.sf(1.0))
  assert_eq(c.sf_at_least(-1.0), 1.0)
  assert_eq(c.quantile(0.01), -1.0)
  assert_eq(c.isf(0.01), 1.0)
  assert_eq(c.quantile(0.5), 0.0)
}

///|
/// 0.3 at the point 0 and 0.7 standard normal: the quantile function has
/// a flat stretch at the atom and follows the normal elsewhere.
test "mixture quantiles by bisection between the components" {
  let m = @dist.Mixture::new([0.3, 0.7], [
    @dist.PointMass::new(0.0) as &@dist.Univariate,
    @dist.Normal::new(),
  ])
  assert_true(rel(m.quantile(0.2), -0.5659488219328631, 1.0e-14))
  assert_eq(m.quantile(0.5), 0.0)
  assert_eq(m.quantile(0.6), 0.0)
  assert_true(rel(m.quantile(0.9), 1.0675705238781412, 1.0e-14))
  assert_true(rel(m.isf(0.1), 1.0675705238781412, 1.0e-14))
  assert_eq(m.quantile(1.0), 1.0 / 0.0)
  let rng = @num.Xoshiro256::seed(3)
  let draws = FixedArray::make(4000, 0.0)
  m.sample_into(rng, draws)
  let mut zeros = 0
  for x in draws {
    if x == 0.0 {
      zeros = zeros + 1
    }
  }
  // 0.3 of 4000, within five standard deviations
  assert_true(zeros > 1055 && zeros < 1345)
}

///|
/// A family holds one parameter fixed and fits the rest.
test "fixing a parameter of a family" {
  let fam = @dist.Exponential::family()
  assert_eq(fam.specs.map(s => s.name), ["loc", "scale"])
  let rate_only = fam.fix("loc", 0.0)
  assert_eq(rate_only.dim(), 1)
  let d = (rate_only.make)([0.5])
  assert_eq((d.loc, d.scale), (0.0, 0.5))
  assert_eq((rate_only.params)(d), [0.5])
  let w = @dist.Weibull::family().fix("loc", 1.0)
  let d = (w.make)([1.5, 2.0])
  assert_eq((d.c, d.loc, d.scale), (1.5, 1.0, 2.0))
  let unknown = try fam.fix("shape", 1.0) catch {
    InvalidParam(name~, ..) => name == "shape"
  } noraise {
    _ => false
  }
  assert_true(unknown)
  let outside = try fam.fix("scale", -1.0) catch {
    InvalidParam(name~, ..) => name == "scale"
  } noraise {
    _ => false
  }
  assert_true(outside)
}

///|
/// A family for a distribution defined outside the package is a record
/// like any other.
test "a family can be written by a caller" {
  let fam : @dist.Family[@dist.PointMass] = {
    name: "point mass",
    specs: [{ name: "at", constraint: Real, }],
    make: p => @dist.PointMass::new(p[0]),
    params: d => [d.at],
  }
  let d = fam.from_unconstrained([2.5])
  assert_eq((fam.params)(d), [2.5])
}
