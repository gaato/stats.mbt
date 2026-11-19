///|
/// Distributions on the integers, written on integer-valued Doubles.
///
/// Each family provides `lattice_pmf`, `lattice_cdf`, `lattice_sf` and
/// `lattice_quantile` (and `lattice_support`), taking the point k as a
/// Double that holds an integer. The generated code in catalog.mbt builds
/// both interfaces on them: the `Double` one (`Cdf`, `Quantile`) passes
/// floor(x) straight through, and the `Int64` one (`Mass`, `IntCdf`,
/// `IntQuantile`) converts at the boundary. No formula here goes through
/// an integer type, so a quantile such as a geometric median of 7e11, or
/// one beyond `Int64`, comes out as it is rather than saturated.

///|
/// 2^53. From here on not every integer is a Double, so a search in steps
/// of one cannot move and a closed-form value is already the answer to the
/// last representable digit.
const EXACT_INTEGERS : Double = 9007199254740992.0

///|
/// The smallest integer k in [lo, hi] with cdf(k) >= q, starting from a
/// guess from a continuous inversion and correcting it by exact
/// comparisons. q = 0 is the lower end of the support and q = 1 the upper.
fn refine_quantile(
  cdf : (Double) -> Double,
  q : Double,
  guess : Double,
  lo : Double,
  hi : Double,
) -> Double {
  if !(q >= 0.0 && q <= 1.0) {
    return 0.0 / 0.0
  }
  if q == 0.0 {
    return lo
  }
  if q == 1.0 {
    return hi
  }
  let g = if guess.is_nan() { lo } else { @math.floor(guess) }
  let mut k = if g < lo { lo } else if g > hi { hi } else { g }
  if k.is_inf() {
    // only a missing guess on an unbounded side gets here; start at 0
    k = if 0.0 < lo { lo } else if 0.0 > hi { hi } else { 0.0 }
  }
  if k.abs() >= EXACT_INTEGERS {
    return k
  }
  while k > lo && cdf(k - 1.0) >= q {
    k = k - 1.0
  }
  while k < hi && cdf(k) < q {
    k = k + 1.0
  }
  k
}

///|
/// The smallest integer k in [lo, hi] with sf(k) <= q, found on the upper
/// tail itself: the bracket is doubled until it holds the answer, then
/// bisected. It needs no guess, and never forms 1 - q.
fn lattice_isf(
  sf : (Double) -> Double,
  q : Double,
  lo : Double,
  hi : Double,
) -> Double {
  if !(q >= 0.0 && q <= 1.0) {
    return 0.0 / 0.0
  }
  if q == 0.0 {
    return hi
  }
  if q == 1.0 {
    return lo
  }
  // invariant: sf(a) > q (or a is below the support), sf(b) <= q
  let mut a = 0.0
  let mut b = 0.0
  let start = if !lo.is_inf() { lo } else if !hi.is_inf() { hi } else { 0.0 }
  if sf(start) <= q {
    if start == lo {
      return lo
    }
    b = start
    let mut step = 1.0
    a = b - step
    while a >= lo && sf(a) <= q {
      b = a
      step = 2.0 * step
      a = b - step
    }
  } else {
    a = start
    let mut step = 1.0
    b = if a + step > hi { hi } else { a + step }
    while b < hi && !b.is_inf() && sf(b) > q {
      a = b
      step = 2.0 * step
      b = if a + step > hi { hi } else { a + step }
    }
  }
  while b - a > 1.0 {
    let m = @math.floor(a + 0.5 * (b - a))
    if m <= a || m >= b {
      break
    }
    if sf(m) <= q {
      b = m
    } else {
      a = m
    }
  }
  b
}

///|
/// An `Int64` from an integer-valued Double, the ends of the range standing
/// for everything beyond them.
fn saturate_int64(x : Double) -> Int64 {
  if x.is_nan() {
    0L
  } else if x >= 9.2233720368547758e18 {
    9223372036854775807L
  } else if x <= -9.2233720368547758e18 {
    -9223372036854775808L
  } else {
    x.to_int64()
  }
}

// ---------------------------------------------------------------------------
// Bernoulli, geometric, discrete uniform
// ---------------------------------------------------------------------------

///|
fn Bernoulli::lattice_pmf(self : Bernoulli, k : Double) -> Double {
  if k == 0.0 {
    1.0 - self.p
  } else if k == 1.0 {
    self.p
  } else {
    0.0
  }
}

///|
fn Bernoulli::lattice_cdf(self : Bernoulli, k : Double) -> Double {
  if k < 0.0 {
    0.0
  } else if k < 1.0 {
    1.0 - self.p
  } else {
    1.0
  }
}

///|
fn Bernoulli::lattice_sf(self : Bernoulli, k : Double) -> Double {
  if k < 0.0 {
    1.0
  } else if k < 1.0 {
    self.p
  } else {
    0.0
  }
}

///|
fn Bernoulli::lattice_quantile(self : Bernoulli, q : Double) -> Double {
  if q.is_nan() {
    q
  } else if q <= 1.0 - self.p {
    0.0
  } else {
    1.0
  }
}

///|
fn Bernoulli::lattice_support(_self : Bernoulli) -> (Double, Double) {
  (0.0, 1.0)
}

///|
fn Geometric::lattice_pmf(self : Geometric, k : Double) -> Double {
  if k < 1.0 {
    0.0
  } else {
    self.p * @math.exp((k - 1.0) * ln1p(-self.p))
  }
}

///|
fn Geometric::lattice_cdf(self : Geometric, k : Double) -> Double {
  if k < 1.0 {
    0.0
  } else {
    -@math.expm1(k * ln1p(-self.p))
  }
}

///|
fn Geometric::lattice_sf(self : Geometric, k : Double) -> Double {
  if k < 1.0 {
    1.0
  } else {
    @math.exp(k * ln1p(-self.p))
  }
}

///|
fn Geometric::lattice_quantile(self : Geometric, q : Double) -> Double {
  let guess = @math.ceil(ln1p(-q) / ln1p(-self.p))
  refine_quantile(k => self.lattice_cdf(k), q, guess, 1.0, 1.0 / 0.0)
}

///|
fn Geometric::lattice_support(_self : Geometric) -> (Double, Double) {
  (1.0, 1.0 / 0.0)
}

///|
fn DiscreteUniform::lattice_pmf(self : DiscreteUniform, k : Double) -> Double {
  let lo = self.low.to_double()
  let hi = self.high.to_double()
  if k < lo || k >= hi {
    0.0
  } else {
    1.0 / (hi - lo)
  }
}

///|
fn DiscreteUniform::lattice_cdf(self : DiscreteUniform, k : Double) -> Double {
  let lo = self.low.to_double()
  let hi = self.high.to_double()
  if k < lo {
    0.0
  } else if k >= hi - 1.0 {
    1.0
  } else {
    (k - lo + 1.0) / (hi - lo)
  }
}

///|
fn DiscreteUniform::lattice_sf(self : DiscreteUniform, k : Double) -> Double {
  let lo = self.low.to_double()
  let hi = self.high.to_double()
  if k < lo {
    1.0
  } else if k >= hi - 1.0 {
    0.0
  } else {
    (hi - 1.0 - k) / (hi - lo)
  }
}

///|
fn DiscreteUniform::lattice_quantile(
  self : DiscreteUniform,
  q : Double,
) -> Double {
  let lo = self.low.to_double()
  let hi = self.high.to_double()
  refine_quantile(
    k => self.lattice_cdf(k),
    q,
    @math.ceil(q * (hi - lo)) + lo - 1.0,
    lo,
    hi - 1.0,
  )
}

///|
fn DiscreteUniform::lattice_support(self : DiscreteUniform) -> (Double, Double) {
  (self.low.to_double(), self.high.to_double() - 1.0)
}

// ---------------------------------------------------------------------------
// Discrete Laplace, Planck, Boltzmann, Zipfian
// ---------------------------------------------------------------------------

///|
fn DiscreteLaplace::lattice_pmf(self : DiscreteLaplace, k : Double) -> Double {
  @math.tanh(0.5 * self.a) * @math.exp(-self.a * k.abs())
}

///|
fn DiscreteLaplace::lattice_cdf(self : DiscreteLaplace, k : Double) -> Double {
  let z = 1.0 + @math.exp(-self.a)
  if k >= 0.0 {
    1.0 - @math.exp(-self.a * (k + 1.0)) / z
  } else {
    @math.exp(self.a * k) / z
  }
}

///|
fn DiscreteLaplace::lattice_sf(self : DiscreteLaplace, k : Double) -> Double {
  let z = 1.0 + @math.exp(-self.a)
  if k >= 0.0 {
    @math.exp(-self.a * (k + 1.0)) / z
  } else {
    1.0 - @math.exp(self.a * k) / z
  }
}

///|
fn DiscreteLaplace::lattice_quantile(
  self : DiscreteLaplace,
  q : Double,
) -> Double {
  let z = 1.0 + @math.exp(-self.a)
  let guess = if q < 1.0 / z {
    @math.floor(@math.ln(q * z) / self.a)
  } else {
    @math.ceil(-@math.ln((1.0 - q) * z) / self.a - 1.0)
  }
  refine_quantile(k => self.lattice_cdf(k), q, guess, -1.0 / 0.0, 1.0 / 0.0)
}

///|
fn DiscreteLaplace::lattice_support(
  _self : DiscreteLaplace,
) -> (Double, Double) {
  (-1.0 / 0.0, 1.0 / 0.0)
}

///|
fn Planck::lattice_pmf(self : Planck, k : Double) -> Double {
  if k < 0.0 {
    0.0
  } else {
    -@math.expm1(-self.lambda) * @math.exp(-self.lambda * k)
  }
}

///|
fn Planck::lattice_cdf(self : Planck, k : Double) -> Double {
  if k < 0.0 {
    0.0
  } else {
    -@math.expm1(-self.lambda * (k + 1.0))
  }
}

///|
fn Planck::lattice_sf(self : Planck, k : Double) -> Double {
  if k < 0.0 {
    1.0
  } else {
    @math.exp(-self.lambda * (k + 1.0))
  }
}

///|
fn Planck::lattice_quantile(self : Planck, q : Double) -> Double {
  refine_quantile(
    k => self.lattice_cdf(k),
    q,
    @math.ceil(-ln1p(-q) / self.lambda - 1.0),
    0.0,
    1.0 / 0.0,
  )
}

///|
fn Planck::lattice_support(_self : Planck) -> (Double, Double) {
  (0.0, 1.0 / 0.0)
}

///|
fn Boltzmann::z(self : Boltzmann) -> Double {
  -@math.expm1(-self.lambda * self.n.to_double())
}

///|
fn Boltzmann::lattice_pmf(self : Boltzmann, k : Double) -> Double {
  if k < 0.0 || k >= self.n.to_double() {
    0.0
  } else {
    -@math.expm1(-self.lambda) * @math.exp(-self.lambda * k) / self.z()
  }
}

///|
fn Boltzmann::lattice_cdf(self : Boltzmann, k : Double) -> Double {
  if k < 0.0 {
    0.0
  } else if k >= self.n.to_double() - 1.0 {
    1.0
  } else {
    -@math.expm1(-self.lambda * (k + 1.0)) / self.z()
  }
}

///|
fn Boltzmann::lattice_sf(self : Boltzmann, k : Double) -> Double {
  let n = self.n.to_double()
  if k < 0.0 {
    1.0
  } else if k >= n - 1.0 {
    0.0
  } else {
    @math.exp(-self.lambda * (k + 1.0)) *
    -@math.expm1(-self.lambda * (n - k - 1.0)) /
    self.z()
  }
}

///|
fn Boltzmann::lattice_quantile(self : Boltzmann, q : Double) -> Double {
  let guess = @math.ceil(-ln1p(-q * self.z()) / self.lambda - 1.0)
  refine_quantile(
    k => self.lattice_cdf(k),
    q,
    guess,
    0.0,
    self.n.to_double() - 1.0,
  )
}

///|
fn Boltzmann::lattice_support(self : Boltzmann) -> (Double, Double) {
  (0.0, self.n.to_double() - 1.0)
}

///|
/// Σ j^(-a) over the integers from..=to.
fn Zipfian::partial(self : Zipfian, from : Double, to : Double) -> Double {
  let mut s = @num.CompensatedSum::zero()
  let mut j = from
  while j <= to {
    s = s.add(@math.pow(j, -self.a))
    j = j + 1.0
  }
  s.value()
}

///|
fn Zipfian::norm(self : Zipfian) -> Double {
  self.partial(1.0, self.n.to_double())
}

///|
fn Zipfian::lattice_pmf(self : Zipfian, k : Double) -> Double {
  if k < 1.0 || k > self.n.to_double() {
    0.0
  } else {
    @math.pow(k, -self.a) / self.norm()
  }
}

///|
fn Zipfian::lattice_cdf(self : Zipfian, k : Double) -> Double {
  if k < 1.0 {
    0.0
  } else if k >= self.n.to_double() {
    1.0
  } else {
    self.partial(1.0, k) / self.norm()
  }
}

///|
/// The upper tail is summed over its own terms, not taken from the lower.
fn Zipfian::lattice_sf(self : Zipfian, k : Double) -> Double {
  let n = self.n.to_double()
  if k < 1.0 {
    1.0
  } else if k >= n {
    0.0
  } else {
    self.partial(k + 1.0, n) / self.norm()
  }
}

///|
fn Zipfian::lattice_quantile(self : Zipfian, q : Double) -> Double {
  refine_quantile(k => self.lattice_cdf(k), q, 1.0, 1.0, self.n.to_double())
}

///|
fn Zipfian::lattice_support(self : Zipfian) -> (Double, Double) {
  (1.0, self.n.to_double())
}
