///|
/// Descriptive statistics as mergeable accumulators. Each is a monoid:
/// `empty()` is the identity and `merge` is associative, so data can be
/// split, summarized in parallel or as a stream, and combined afterwards.
/// Values are immutable so the types can later become `#valtype` on native.
pub(open) trait Accumulator {
  empty() -> Self
  push(Self, Double) -> Self
  merge(Self, Self) -> Self
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
pub impl Accumulator for MeanVar with empty() {
  { n: 0, mean: 0.0, m2: 0.0 }
}

///|
pub impl Accumulator for MeanVar with push(self, x) {
  let n = self.n + 1
  let d = x - self.mean
  let mean = self.mean + d / n.to_double()
  { n, mean, m2: self.m2 + d * (x - mean) }
}

///|
pub impl Accumulator for MeanVar with merge(a, b) {
  if a.n == 0 {
    return b
  }
  if b.n == 0 {
    return a
  }
  let n = a.n + b.n
  let d = b.mean - a.mean
  let fb = b.n.to_double() / n.to_double()
  {
    n,
    mean: a.mean + d * fb,
    m2: a.m2 + b.m2 + d * d * a.n.to_double() * fb,
  }
}

///|
pub fn MeanVar::count(self : MeanVar) -> Int64 {
  self.n
}

///|
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
pub fn MeanVar::variance(self : MeanVar, ddof~ : Int) -> Double? {
  let d = self.n - ddof.to_int64()
  if d <= 0 {
    None
  } else {
    Some(self.m2 / d.to_double())
  }
}

///|
/// Trait methods become ordinary methods for callers in other packages.
pub extend MeanVar with Accumulator::{push, merge}
