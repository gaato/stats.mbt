///|
/// The two test cases for the capability traits besides the parametric
/// families: a distribution built from data, and one built from a table of
/// probabilities (for example an exact rank-test null distribution obtained
/// by polynomial multiplication or FFT convolution).

///|
pub struct Empirical {
  priv sorted : FixedArray[Double]
  priv quantile_method : @data.QuantileMethod
}

///|
/// Rejects empty data and NaN instead of dropping them silently.
pub fn Empirical::new(
  data : ArrayView[Double],
  quantile_method~ : @data.QuantileMethod,
) -> Empirical raise ParamError {
  if data.length() == 0 {
    raise InvalidParam(name="data", value=0.0, reason="must not be empty")
  }
  for x in data {
    if x.is_nan() {
      raise InvalidParam(name="data", value=x, reason="contains NaN")
    }
  }
  let sorted = FixedArray::from_iter(data.iter())
  sorted.sort()
  { sorted, quantile_method }
}

///|
/// The number of values <= x, by binary search.
fn Empirical::count_le(self : Empirical, x : Double) -> Int {
  let mut lo = 0
  let mut hi = self.sorted.length()
  while lo < hi {
    let mid = (lo + hi) / 2
    if self.sorted[mid] <= x {
      lo = mid + 1
    } else {
      hi = mid
    }
  }
  lo
}

///|
pub impl Cdf for Empirical with cdf(self, x) {
  self.count_le(x).to_double() / self.sorted.length().to_double()
}

///|
pub impl Cdf for Empirical with sf(self, x) {
  (self.sorted.length() - self.count_le(x)).to_double() /
  self.sorted.length().to_double()
}

///|
/// P(X >= x): the values not below x.
pub impl Cdf for Empirical with sf_at_least(self, x) {
  let mut lo = 0
  let mut hi = self.sorted.length()
  while lo < hi {
    let mid = (lo + hi) / 2
    if self.sorted[mid] < x {
      lo = mid + 1
    } else {
      hi = mid
    }
  }
  (self.sorted.length() - lo).to_double() / self.sorted.length().to_double()
}

///|
pub impl Quantile for Empirical with quantile(self, p) {
  @describe.quantile_sorted(self.sorted[:], p, self.quantile_method)
}

///|
/// Bootstrap resampling is sampling from the empirical distribution.
pub impl Sample for Empirical with sample(self, rng) {
  self.sorted[@num.uniform_below(rng, self.sorted.length().to_uint64()).to_int()]
}

///|
pub impl Moments for Empirical with mean(self) {
  let mut s = @num.CompensatedSum::zero()
  for x in self.sorted {
    s = s.add(x)
  }
  Some(s.value() / self.sorted.length().to_double())
}

///|
/// The variance of the distribution itself: divisor n.
pub impl Moments for Empirical with variance(self) {
  let m = Moments::mean(self).unwrap()
  let mut s = @num.CompensatedSum::zero()
  for x in self.sorted {
    s = s.add((x - m) * (x - m))
  }
  Some(s.value() / self.sorted.length().to_double())
}

///|
/// Probabilities on the integers `offset ..< offset + probs.length()`.
/// Both prefix and suffix sums are kept, so `sf` far in the upper tail
/// is a sum of small terms, not 1 minus a number close to 1.
pub struct FiniteDiscrete {
  priv offset : Int
  priv probs : FixedArray[Double]
  priv lower : FixedArray[Double]
  priv upper : FixedArray[Double]
}

///|
pub fn FiniteDiscrete::new(
  offset : Int,
  probs : ArrayView[Double],
) -> FiniteDiscrete raise ParamError {
  let n = probs.length()
  if n == 0 {
    raise InvalidParam(name="probs", value=0.0, reason="must not be empty")
  }
  let mut total = @num.CompensatedSum::zero()
  for q in probs {
    if !(q >= 0.0) {
      raise InvalidParam(name="probs", value=q, reason="must be >= 0")
    }
    total = total.add(q)
  }
  if (total.value() - 1.0).abs() > 1.0e-9 {
    raise InvalidParam(name="probs", value=total.value(), reason="must sum to 1")
  }
  let lower = FixedArray::make(n, 0.0)
  let upper = FixedArray::make(n, 0.0)
  let mut acc = @num.CompensatedSum::zero()
  for i in 0..<n {
    acc = acc.add(probs[i])
    lower[i] = acc.value()
  }
  acc = @num.CompensatedSum::zero()
  for i = n - 1; i >= 0; i = i - 1 {
    // upper[i] = P(X > offset + i)
    upper[i] = acc.value()
    acc = acc.add(probs[i])
  }
  { offset, probs: FixedArray::from_iter(probs.iter()), lower, upper }
}

///|
pub impl Mass for FiniteDiscrete with pmf(self, k) {
  let i = k - self.offset
  if i < 0 || i >= self.probs.length() {
    0.0
  } else {
    self.probs[i]
  }
}

///|
pub impl IntCdf for FiniteDiscrete with cdf_int(self, k) {
  let i = k - self.offset
  if i < 0 {
    0.0
  } else if i >= self.probs.length() {
    1.0
  } else {
    self.lower[i]
  }
}

///|
pub impl IntCdf for FiniteDiscrete with sf_int(self, k) {
  let i = k - self.offset
  if i < 0 {
    1.0
  } else if i >= self.probs.length() {
    0.0
  } else {
    self.upper[i]
  }
}

///|
pub impl IntQuantile for FiniteDiscrete with quantile_int(self, p) {
  let mut lo = 0
  let mut hi = self.lower.length() - 1
  while lo < hi {
    let mid = (lo + hi) / 2
    if self.lower[mid] >= p {
      hi = mid
    } else {
      lo = mid + 1
    }
  }
  self.offset + lo
}

///|
pub impl Quantile for FiniteDiscrete with quantile(self, p) {
  IntQuantile::quantile_int(self, p).to_double()
}

///|
pub impl Sample for FiniteDiscrete with sample(self, rng) {
  Quantile::quantile(self, @num.uniform_open(rng))
}

///|
pub impl Univariate for FiniteDiscrete

///|
pub impl Univariate for Empirical

///|
/// On the real line the distribution function is the step function.
pub impl Cdf for FiniteDiscrete with cdf(self, x) {
  IntCdf::cdf_int(self, @math.floor(x).to_int())
}

///|
pub impl Cdf for FiniteDiscrete with sf(self, x) {
  IntCdf::sf_int(self, @math.floor(x).to_int())
}

///|
pub impl Cdf for FiniteDiscrete with sf_at_least(self, x) {
  IntCdf::sf_int(self, @math.ceil(x).to_int() - 1)
}
