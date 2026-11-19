///|
/// Descriptive statistics as mergeable accumulators. Each is a monoid:
/// `empty()` is the identity and `merge` is associative, so data can be
/// split, summarized in parallel or as a stream, and combined afterwards.
/// Values are immutable so the types can later become `#valtype` on native.
pub(open) trait Accumulator {
  fn empty() -> Self
  fn push(Self, Double) -> Self
  fn merge(Self, Self) -> Self
}

///|
/// Count, mean and the sum of squared deviations (Welford). `merge` uses
/// the pairwise update of Chan, Golub and LeVeque, which stays stable when
/// the two halves have very different means.
pub struct MeanVar {
  priv n : Int64
  priv mean : Double
  priv m2 : Double
}

///|
/// The accumulator of no values.
pub impl Accumulator for MeanVar with fn empty() {
  { n: 0, mean: 0.0, m2: 0.0, }
}

///|
/// Adds one value (Welford's update).
pub impl Accumulator for MeanVar with fn push(self, x) {
  let n = self.n + 1
  let d = x - self.mean
  let mean = self.mean + d / n.to_double()
  { n, mean, m2: self.m2 + d * (x - mean), }
}

///|
/// Combines two accumulators as if their values had been pushed into one
/// (Chan et al.).
pub impl Accumulator for MeanVar with fn merge(a, b) {
  if a.n == 0 {
    return b
  }
  if b.n == 0 {
    return a
  }
  let n = a.n + b.n
  let d = b.mean - a.mean
  let fb = b.n.to_double() / n.to_double()
  { n, mean: a.mean + d * fb, m2: a.m2 + b.m2 + d * d * a.n.to_double() * fb, }
}

///|
/// The number of values pushed.
pub fn MeanVar::count(self : MeanVar) -> Int64 {
  self.n
}

///|
/// The mean, or `None` before any value.
pub fn MeanVar::mean(self : MeanVar) -> Double? {
  if self.n == 0 {
    None
  } else {
    Some(self.mean)
  }
}

///|
/// The divisor is an argument, not a default: ddof=1 is the sample
/// variance (R, pandas), ddof=0 the population variance (NumPy's default).
pub fn MeanVar::variance(self : MeanVar, ddof~ : Int32) -> Double? {
  let d = self.n - ddof.to_int64()
  if d <= 0 {
    None
  } else {
    Some(self.m2 / d.to_double())
  }
}

///|
/// Trait methods become ordinary methods for callers in other packages:
/// `MeanVar::empty().push(x)`.
pub extend MeanVar with Accumulator::{empty, push, merge}

///|
/// Count, smallest and largest value. NaN is skipped by the comparisons,
/// so it is counted but never becomes an extreme; clean the data first if
/// that is not wanted.
pub struct Extremes {
  priv n : Int64
  priv lo : Double
  priv hi : Double
}

///|
/// The accumulator of no values.
pub impl Accumulator for Extremes with fn empty() {
  { n: 0, lo: 1.0 / 0.0, hi: -1.0 / 0.0, }
}

///|
/// Adds one value.
pub impl Accumulator for Extremes with fn push(self, x) {
  {
    n: self.n + 1,
    lo: if x < self.lo {
      x
    } else {
      self.lo
    },
    hi: if x > self.hi {
      x
    } else {
      self.hi
    },
  }
}

///|
/// Combines two accumulators as if their values had been pushed into one.
pub impl Accumulator for Extremes with fn merge(a, b) {
  {
    n: a.n + b.n,
    lo: if b.lo < a.lo {
      b.lo
    } else {
      a.lo
    },
    hi: if b.hi > a.hi {
      b.hi
    } else {
      a.hi
    },
  }
}

///|
/// Dot syntax for `Accumulator`: `d.empty`, `d.push`, `d.merge`.
pub extend Extremes with Accumulator::{empty, push, merge}

///|
/// The number of values pushed.
pub fn Extremes::count(self : Extremes) -> Int64 {
  self.n
}

///|
/// `None` before any value has been pushed.
pub fn Extremes::min(self : Extremes) -> Double? {
  if self.n == 0 {
    None
  } else {
    Some(self.lo)
  }
}

///|
/// The largest value, or `None` before any value.
pub fn Extremes::max(self : Extremes) -> Double? {
  if self.n == 0 {
    None
  } else {
    Some(self.hi)
  }
}

///|
/// One pass over an array with any accumulator.
pub fn[A : Accumulator] accumulate(data : ArrayView[Double]) -> A {
  let mut acc : A = Accumulator::empty()
  for x in data {
    acc = acc.push(x)
  }
  acc
}
