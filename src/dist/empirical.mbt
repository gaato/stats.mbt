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
pub impl Cdf for Empirical with cdf(self, x) {
  ignore((self, x))
  ...
}

///|
pub impl Quantile for Empirical with quantile(self, p) {
  ignore((self, p))
  ...
}

///|
/// Bootstrap resampling is sampling from the empirical distribution.
pub impl Sample for Empirical with sample(self, rng) {
  ignore((self, rng))
  ...
}

///|
pub impl Moments for Empirical with mean(self) {
  ignore(self)
  ...
}

///|
pub impl Moments for Empirical with variance(self) {
  ignore(self)
  ...
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
  ignore((offset, probs))
  ...
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
  ignore((self, p))
  ...
}

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
