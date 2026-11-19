///|
/// Distributions that are not parametric families: one built from data,
/// one from a table of probabilities on the integers (for example an exact
/// rank-test null distribution obtained by polynomial multiplication or
/// FFT convolution), and one from atoms anywhere on the real line.

///|
/// The distribution that puts mass 1/n on each of n observed values.
#derive(Debug)
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
  { sorted, quantile_method, }
}

///|
/// The observed values, ascending: what a histogram of a bootstrap or
/// permutation distribution is drawn from.
pub fn Empirical::values(self : Empirical) -> ArrayView[Double] {
  self.sorted[:]
}

///|
/// The number of observed values.
pub fn Empirical::length(self : Empirical) -> Int32 {
  self.sorted.length()
}

///|
/// The number of entries of an ascending array that are <= x (`strict`
/// false) or < x (`strict` true), by binary search.
fn count_below(
  sorted : FixedArray[Double],
  x : Double,
  strict~ : Bool,
) -> Int32 {
  let mut lo = 0
  let mut hi = sorted.length()
  while lo < hi {
    let mid = (lo + hi) / 2
    if sorted[mid] < x || (!strict && sorted[mid] == x) {
      lo = mid + 1
    } else {
      hi = mid
    }
  }
  lo
}

///|
/// P(X = x): the share of observed values equal to x.
pub fn Empirical::mass(self : Empirical, x : Double) -> Double {
  let n = self.sorted.length()
  (count_below(self.sorted, x, strict=false) -
  count_below(self.sorted, x, strict=true)).to_double() /
  n.to_double()
}

///|
/// P(X ≤ x).
pub impl Cdf for Empirical with fn cdf(self, x) {
  count_below(self.sorted, x, strict=false).to_double() /
  self.sorted.length().to_double()
}

///|
/// Counted from the top, not taken from the lower tail.
pub impl Cdf for Empirical with fn sf(self, x) {
  let n = self.sorted.length()
  (n - count_below(self.sorted, x, strict=false)).to_double() / n.to_double()
}

///|
/// P(X >= x): the values not below x.
pub impl Cdf for Empirical with fn sf_at_least(self, x) {
  let n = self.sorted.length()
  (n - count_below(self.sorted, x, strict=true)).to_double() / n.to_double()
}

///|
/// The sample quantile by the definition chosen at construction.
pub impl Quantile for Empirical with fn quantile(self, p) {
  @describe.quantile_sorted(self.sorted[:], p, self.quantile_method)
}

///|
/// The sample quantile at 1 - q. Unlike for a parametric family this
/// loses nothing: a sample of n values has no quantile finer than 1/n, and
/// below that every definition returns the largest value.
pub impl Quantile for Empirical with fn isf(self, q) {
  @describe.quantile_sorted(self.sorted[:], 1.0 - q, self.quantile_method)
}

///|
/// Bootstrap resampling is sampling from the empirical distribution.
pub impl Sample for Empirical with fn sample(self, rng) {
  self.sorted[@num.uniform_index(rng, self.sorted.length())]
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl Univariate for Empirical

///|
/// The mean, or `None` when it does not exist.
pub impl Moments for Empirical with fn mean(self) {
  let mut s = @num.CompensatedSum::zero()
  for x in self.sorted {
    s = s.add(x)
  }
  Some(s.value() / self.sorted.length().to_double())
}

///|
/// The variance of the distribution itself: divisor n.
pub impl Moments for Empirical with fn variance(self) {
  let m = Moments::mean(self).unwrap()
  let mut s = @num.CompensatedSum::zero()
  for x in self.sorted {
    s = s.add((x - m) * (x - m))
  }
  Some(s.value() / self.sorted.length().to_double())
}

///|
/// The lower end of the support.
pub impl Support for Empirical with fn support_lower(self) {
  self.sorted[0]
}

///|
/// The upper end of the support.
pub impl Support for Empirical with fn support_upper(self) {
  self.sorted[self.sorted.length() - 1]
}

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend Empirical with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend Empirical with Quantile::{quantile, isf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend Empirical with Sample::{sample, sample_into}

///|
/// Dot syntax for `Moments`: `d.mean`, `d.variance`.
pub extend Empirical with Moments::{mean, variance}

///|
/// Dot syntax for `Support`: `d.support_lower`, `d.support_upper`.
pub extend Empirical with Support::{support_lower, support_upper}

///|
/// Running sums of probabilities from both ends: `lower[i]` is the total
/// up to and including entry i, `upper[i]` the total after it. Each tail
/// is then a sum of its own terms, not 1 minus the other.
fn running_sums(
  probs : FixedArray[Double],
) -> (FixedArray[Double], FixedArray[Double]) {
  let n = probs.length()
  let lower = FixedArray::make(n, 0.0)
  let upper = FixedArray::make(n, 0.0)
  let mut acc = @num.CompensatedSum::zero()
  for i in 0..<n {
    acc = acc.add(probs[i])
    lower[i] = acc.value()
  }
  acc = @num.CompensatedSum::zero()
  for i = n - 1; i >= 0; i = i - 1 {
    upper[i] = acc.value()
    acc = acc.add(probs[i])
  }
  (lower, upper)
}

///|
/// The first index where `reached` holds, for a condition that is false up
/// to some entry and true from it on; the last index if it never holds.
fn first_index(n : Int32, reached : (Int32) -> Bool) -> Int32 {
  let mut lo = 0
  let mut hi = n - 1
  while lo < hi {
    let mid = (lo + hi) / 2
    if reached(mid) {
      hi = mid
    } else {
      lo = mid + 1
    }
  }
  lo
}

///|
/// Probabilities on the integers `offset ..< offset + probs.length()`.
/// Both prefix and suffix sums are kept, so `sf` far in the upper tail
/// is a sum of small terms, not 1 minus a number close to 1.
#derive(Debug)
pub struct FiniteDiscrete {
  priv offset : Int64
  priv probs : FixedArray[Double]
  priv lower : FixedArray[Double]
  priv upper : FixedArray[Double]
}

///|
/// Probabilities of the integers offset, offset + 1, …; they must be
/// non-negative and sum to 1.
pub fn FiniteDiscrete::new(
  offset : Int64,
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
    raise InvalidParam(
      name="probs",
      value=total.value(),
      reason="must sum to 1",
    )
  }
  let probs = FixedArray::from_iter(probs.iter())
  let (lower, upper) = running_sums(probs)
  { offset, probs, lower, upper, }
}

///|
/// The integer the first probability belongs to.
pub fn FiniteDiscrete::offset(self : FiniteDiscrete) -> Int64 {
  self.offset
}

///|
/// The probabilities of `offset`, `offset + 1`, and so on.
pub fn FiniteDiscrete::probs(self : FiniteDiscrete) -> ArrayView[Double] {
  self.probs[:]
}

///|
/// The position of the integer-valued Double k in the table, as a Double
/// so that a point far outside it needs no integer conversion.
fn FiniteDiscrete::position(self : FiniteDiscrete, k : Double) -> Double {
  k - self.offset.to_double()
}

///|
fn FiniteDiscrete::lattice_pmf(self : FiniteDiscrete, k : Double) -> Double {
  let i = self.position(k)
  if i < 0.0 || i >= self.probs.length().to_double() {
    0.0
  } else {
    self.probs[i.to_int()]
  }
}

///|
fn FiniteDiscrete::lattice_cdf(self : FiniteDiscrete, k : Double) -> Double {
  let i = self.position(k)
  if i < 0.0 {
    0.0
  } else if i >= self.probs.length().to_double() {
    1.0
  } else {
    self.lower[i.to_int()]
  }
}

///|
fn FiniteDiscrete::lattice_sf(self : FiniteDiscrete, k : Double) -> Double {
  let i = self.position(k)
  if i < 0.0 {
    1.0
  } else if i >= self.probs.length().to_double() {
    0.0
  } else {
    self.upper[i.to_int()]
  }
}

///|
/// P(X = k).
pub impl Mass for FiniteDiscrete with fn pmf(self, k) {
  self.lattice_pmf(k.to_double())
}

///|
/// P(X ≤ k).
pub impl IntCdf for FiniteDiscrete with fn cdf_int(self, k) {
  self.lattice_cdf(k.to_double())
}

///|
/// P(X > k).
pub impl IntCdf for FiniteDiscrete with fn sf_int(self, k) {
  self.lattice_sf(k.to_double())
}

///|
/// The smallest k with P(X ≤ k) ≥ p.
pub impl IntQuantile for FiniteDiscrete with fn quantile_int(self, p) {
  self.offset +
  first_index(self.lower.length(), i => self.lower[i] >= p).to_int64()
}

///|
/// One draw from `rng`.
pub impl SampleInt for FiniteDiscrete with fn sample_int(self, rng) {
  IntQuantile::quantile_int(self, @num.uniform_open(rng))
}

///|
/// On the real line the distribution function is the step function.
pub impl Cdf for FiniteDiscrete with fn cdf(self, x) {
  self.lattice_cdf(@math.floor(x))
}

///|
/// P(X > x).
pub impl Cdf for FiniteDiscrete with fn sf(self, x) {
  self.lattice_sf(@math.floor(x))
}

///|
/// P(X ≥ x).
pub impl Cdf for FiniteDiscrete with fn sf_at_least(self, x) {
  self.lattice_sf(@math.ceil(x) - 1.0)
}

///|
/// The smallest x with P(X ≤ x) ≥ p.
pub impl Quantile for FiniteDiscrete with fn quantile(self, p) {
  if p.is_nan() {
    return p
  }
  self.offset.to_double() +
  first_index(self.lower.length(), i => self.lower[i] >= p).to_double()
}

///|
/// The smallest k with P(X > k) <= q, on the suffix sums.
pub impl Quantile for FiniteDiscrete with fn isf(self, q) {
  if q.is_nan() {
    return q
  }
  self.offset.to_double() +
  first_index(self.upper.length(), i => self.upper[i] <= q).to_double()
}

///|
/// One draw from `rng`.
pub impl Sample for FiniteDiscrete with fn sample(self, rng) {
  Quantile::quantile(self, @num.uniform_open(rng))
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl Univariate for FiniteDiscrete

///|
/// The first integer with positive probability.
pub impl Support for FiniteDiscrete with fn support_lower(self) {
  let mut i = 0
  while i < self.probs.length() - 1 && self.probs[i] == 0.0 {
    i = i + 1
  }
  self.offset.to_double() + i.to_double()
}

///|
/// The last integer with positive probability.
pub impl Support for FiniteDiscrete with fn support_upper(self) {
  let mut i = self.probs.length() - 1
  while i > 0 && self.probs[i] == 0.0 {
    i = i - 1
  }
  self.offset.to_double() + i.to_double()
}

///|
/// Dot syntax for `Mass`: `d.pmf`, `d.logpmf`.
pub extend FiniteDiscrete with Mass::{pmf, logpmf}

///|
/// Dot syntax for `IntCdf`: `d.cdf_int`, `d.sf_int`.
pub extend FiniteDiscrete with IntCdf::{cdf_int, sf_int}

///|
/// Dot syntax for `IntQuantile`: `d.quantile_int`.
pub extend FiniteDiscrete with IntQuantile::{quantile_int}

///|
/// Dot syntax for `SampleInt`: `d.sample_int`.
pub extend FiniteDiscrete with SampleInt::{sample_int}

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend FiniteDiscrete with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend FiniteDiscrete with Quantile::{quantile, isf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend FiniteDiscrete with Sample::{sample, sample_into}

///|
/// Dot syntax for `Support`: `d.support_lower`, `d.support_upper`.
pub extend FiniteDiscrete with Support::{support_lower, support_upper}

///|
/// Finitely many atoms anywhere on the real line, each with its
/// probability: the exact null distribution of a rank statistic with ties
/// (half-integers), a weighted sample, or an `Empirical` with its repeated
/// values merged. Unlike `Empirical` it can say how much mass sits on a
/// point, which the minimum-likelihood two-sided p-value needs.
///
/// Values are kept ascending and distinct, with prefix and suffix sums, so
/// both tails are exact sums of the atoms on their own side.
#derive(Debug)
pub struct Atoms {
  priv values : FixedArray[Double]
  priv probs : FixedArray[Double]
  priv lower : FixedArray[Double]
  priv upper : FixedArray[Double]
}

///|
/// From value and probability pairs that are already valid: sorts, merges
/// equal values, drops atoms without mass and builds the running sums.
fn Atoms::build(pairs : Array[(Double, Double)]) -> Atoms {
  pairs.sort_by((a, b) => a.0.compare(b.0))
  let values = []
  let probs = []
  for pair in pairs {
    let (v, q) = pair
    if q == 0.0 {
      continue
    }
    let n = values.length()
    if n > 0 && values[n - 1] == v {
      probs[n - 1] = probs[n - 1] + q
    } else {
      values.push(v)
      probs.push(q)
    }
  }
  let probs = FixedArray::from_iter(probs.iter())
  let (lower, upper) = running_sums(probs)
  { values: FixedArray::from_iter(values.iter()), probs, lower, upper, }
}

///|
/// Atoms at `values` with probabilities `probs`, in any order; a value
/// given twice has its probabilities added. The probabilities must be
/// non-negative and sum to 1 within rounding.
pub fn Atoms::new(
  values : ArrayView[Double],
  probs : ArrayView[Double],
) -> Atoms raise ParamError {
  if values.length() != probs.length() || values.length() == 0 {
    raise InvalidParam(
      name="probs",
      value=probs.length().to_double(),
      reason="need one probability per value, at least one",
    )
  }
  let mut total = @num.CompensatedSum::zero()
  let pairs = []
  for i, v in values {
    if v.is_nan() || v.is_inf() {
      raise InvalidParam(name="values", value=v, reason="must be finite")
    }
    if !(probs[i] >= 0.0) {
      raise InvalidParam(name="probs", value=probs[i], reason="must be >= 0")
    }
    total = total.add(probs[i])
    pairs.push((v, probs[i]))
  }
  if (total.value() - 1.0).abs() > 1.0e-9 {
    raise InvalidParam(
      name="probs",
      value=total.value(),
      reason="must sum to 1",
    )
  }
  Atoms::build(pairs)
}

///|
/// The distribution of a weighted sample: each value with its weight over
/// the total weight. The name says that the weights are normalised; `new`
/// never rescales what it is given.
pub fn Atoms::from_weighted(w : @data.Weighted) -> Atoms raise ParamError {
  let mut total = @num.CompensatedSum::zero()
  for i, v in w.values {
    if v.is_nan() || v.is_inf() {
      raise InvalidParam(name="values", value=v, reason="must be finite")
    }
    total = total.add(w.weights.values[i])
  }
  let sum = total.value()
  if !(sum > 0.0) || sum.is_inf() {
    raise InvalidParam(
      name="weights",
      value=sum,
      reason="must have a positive, finite total",
    )
  }
  let pairs = []
  for i, v in w.values {
    pairs.push((v, w.weights.values[i] / sum))
  }
  Atoms::build(pairs)
}

///|
/// The same distribution with repeated values merged into one atom each.
pub fn Empirical::to_atoms(self : Empirical) -> Atoms {
  let share = 1.0 / self.sorted.length().to_double()
  let pairs = []
  let mut i = 0
  while i < self.sorted.length() {
    let mut j = i
    while j + 1 < self.sorted.length() && self.sorted[j + 1] == self.sorted[i] {
      j = j + 1
    }
    pairs.push((self.sorted[i], (j - i + 1).to_double() * share))
    i = j + 1
  }
  Atoms::build(pairs)
}

///|
/// The atoms, ascending and distinct.
pub fn Atoms::values(self : Atoms) -> ArrayView[Double] {
  self.values[:]
}

///|
/// The probability of each atom, in the order of `values`.
pub fn Atoms::probs(self : Atoms) -> ArrayView[Double] {
  self.probs[:]
}

///|
/// The number of atoms.
pub fn Atoms::length(self : Atoms) -> Int32 {
  self.values.length()
}

///|
/// P(X = x): the probability of the atom at x, or 0.
pub fn Atoms::mass(self : Atoms, x : Double) -> Double {
  let i = count_below(self.values, x, strict=true)
  if i < self.values.length() && self.values[i] == x {
    self.probs[i]
  } else {
    0.0
  }
}

///|
/// P(X ≤ x).
pub impl Cdf for Atoms with fn cdf(self, x) {
  let k = count_below(self.values, x, strict=false)
  if k == 0 {
    0.0
  } else if k == self.values.length() {
    1.0
  } else {
    self.lower[k - 1]
  }
}

///|
/// P(X > x).
pub impl Cdf for Atoms with fn sf(self, x) {
  let k = count_below(self.values, x, strict=false)
  if k == 0 {
    1.0
  } else {
    self.upper[k - 1]
  }
}

///|
/// P(X >= x): everything from the first atom not below x.
pub impl Cdf for Atoms with fn sf_at_least(self, x) {
  let k = count_below(self.values, x, strict=true)
  if k == 0 {
    1.0
  } else {
    self.upper[k - 1]
  }
}

///|
/// The smallest x with P(X ≤ x) ≥ p.
pub impl Quantile for Atoms with fn quantile(self, p) {
  if p.is_nan() {
    return p
  }
  self.values[first_index(self.lower.length(), i => self.lower[i] >= p)]
}

///|
/// The smallest x with P(X > x) ≤ q.
pub impl Quantile for Atoms with fn isf(self, q) {
  if q.is_nan() {
    return q
  }
  self.values[first_index(self.upper.length(), i => self.upper[i] <= q)]
}

///|
/// One draw from `rng`.
pub impl Sample for Atoms with fn sample(self, rng) {
  Quantile::quantile(self, @num.uniform_open(rng))
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl Univariate for Atoms

///|
/// The mean, or `None` when it does not exist.
pub impl Moments for Atoms with fn mean(self) {
  let mut s = @num.CompensatedSum::zero()
  for i, v in self.values {
    s = s.add(self.probs[i] * v)
  }
  Some(s.value())
}

///|
/// The variance, or `None` when it does not exist.
pub impl Moments for Atoms with fn variance(self) {
  let m = Moments::mean(self).unwrap()
  let mut s = @num.CompensatedSum::zero()
  for i, v in self.values {
    s = s.add(self.probs[i] * (v - m) * (v - m))
  }
  Some(s.value())
}

///|
/// The lower end of the support.
pub impl Support for Atoms with fn support_lower(self) {
  self.values[0]
}

///|
/// The upper end of the support.
pub impl Support for Atoms with fn support_upper(self) {
  self.values[self.values.length() - 1]
}

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend Atoms with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend Atoms with Quantile::{quantile, isf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend Atoms with Sample::{sample, sample_into}

///|
/// Dot syntax for `Moments`: `d.mean`, `d.variance`.
pub extend Atoms with Moments::{mean, variance}

///|
/// Dot syntax for `Support`: `d.support_lower`, `d.support_upper`.
pub extend Atoms with Support::{support_lower, support_upper}
