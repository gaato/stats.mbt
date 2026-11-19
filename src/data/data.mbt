///|
/// Data shapes below the statistics: what one observation is, how missing
/// values are treated, what weights mean, and how counts are tabulated.
/// Nothing here knows about distributions, so every other package can
/// depend on it.

///|
/// One observation as it enters a likelihood. Censoring is a property of
/// how a value was observed, not of the distribution: a right-censored
/// survival time contributes log S(x), an interval-censored one
/// log(F(hi) - F(lo)). (A `Censored` distribution in `dist` is a different
/// thing: values clamped at detection limits.)
#derive(Eq, Debug)
pub(all) enum Observation {
  /// The value itself.
  Exact(Double)
  /// The value is greater than this.
  RightCensored(Double)
  /// The value is at most this.
  LeftCensored(Double)
  /// lo < value <= hi.
  IntervalCensored(lo~ : Double, hi~ : Double)
}

///|
/// What to do with NaN in the input, chosen by the caller every time.
/// `Reject` raises; `Omit` drops and reports how many were dropped;
/// `Propagate` returns NaN results. (SciPy's nan_policy, R's na.rm.)
#derive(Eq, Debug)
pub(all) enum MissingPolicy {
  Reject
  Omit
  Propagate
}

///|
#derive(Debug)
pub(all) suberror DataError {
  ContainsNaN(index~ : Int)
  LengthMismatch(expected~ : Int, got~ : Int)
  Empty
  NegativeWeight(index~ : Int)
  NegativeCount(row~ : Int, col~ : Int)
}

///|
/// The meaning of weights changes the variance, so it is part of the data:
/// frequency weights count repeated rows (n = Σw); reliability weights are
/// inverse variances (the effective n is (Σw)²/Σw²); sampling weights are
/// inverse inclusion probabilities from a survey design.
#derive(Eq, Debug)
pub(all) enum WeightKind {
  Frequency
  Reliability
  Sampling
}

///|
pub struct Weighted {
  values : FixedArray[Double]
  weights : FixedArray[Double]
  kind : WeightKind
}

///|
pub fn Weighted::new(
  values : ArrayView[Double],
  weights : ArrayView[Double],
  kind~ : WeightKind,
) -> Weighted raise DataError {
  if values.length() != weights.length() {
    raise LengthMismatch(expected=values.length(), got=weights.length())
  }
  for i, w in weights {
    if !(w >= 0.0) {
      raise NegativeWeight(index=i)
    }
  }
  {
    values: FixedArray::from_iter(values.iter()),
    weights: FixedArray::from_iter(weights.iter()),
    kind,
  }
}

///|
/// An r × c table of counts, row-major. Tests of independence, Fisher's
/// exact test and measures of association take this, not raw labels.
pub struct Table {
  rows : Int
  cols : Int
  counts : FixedArray[Int64]
}

///|
pub fn Table::new(
  rows : Int,
  cols : Int,
  counts : ArrayView[Int64],
) -> Table raise DataError {
  if counts.length() != rows * cols {
    raise LengthMismatch(expected=rows * cols, got=counts.length())
  }
  for i, c in counts {
    if c < 0 {
      raise NegativeCount(row=i / cols, col=i % cols)
    }
  }
  { rows, cols, counts: FixedArray::from_iter(counts.iter()) }
}

///|
pub fn Table::at(self : Table, row : Int, col : Int) -> Int64 {
  self.counts[row * self.cols + col]
}

///|
pub fn Table::total(self : Table) -> Int64 {
  let mut s = 0L
  for c in self.counts {
    s = s + c
  }
  s
}

///|
/// Apply a missing-value policy to one sample. Returns the cleaned values
/// and the number omitted; `Propagate` keeps NaN in place.
pub fn clean(
  values : ArrayView[Double],
  policy : MissingPolicy,
) -> (FixedArray[Double], Int) raise DataError {
  match policy {
    Propagate => (FixedArray::from_iter(values.iter()), 0)
    Reject => {
      for i, x in values {
        if x.is_nan() {
          raise ContainsNaN(index=i)
        }
      }
      (FixedArray::from_iter(values.iter()), 0)
    }
    Omit => {
      let kept = values.iter().filter(x => !x.is_nan()).collect()
      (FixedArray::from_iter(kept.iter()), values.length() - kept.length())
    }
  }
}
