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
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum MissingPolicy {
  Reject
  Omit
  Propagate
}

///|
/// Why data could not be accepted: what was wrong and where.
#non_exhaustive
#derive(Debug)
pub(all) suberror DataError {
  ContainsNaN(index~ : Int32)
  LengthMismatch(expected~ : Int32, got~ : Int32)
  Empty
  /// A weight that is negative or NaN.
  NegativeWeight(index~ : Int32)
  NegativeCount(row~ : Int32, col~ : Int32)
  /// Row and column counts must both be positive.
  InvalidShape(rows~ : Int32, cols~ : Int32)
  /// A group, row or column label outside `0 ..< count`.
  InvalidLabel(index~ : Int32, label~ : Int32)
  /// A probability or proportion outside the range the function accepts.
  InvalidProbability(Double)
  /// Fewer values than the statistic is defined for (a variance with
  /// divisor n - 1 needs two).
  TooFewValues(needed~ : Int32, got~ : Int32)
}

///|
/// The meaning of weights changes the variance, so it is part of the data:
/// frequency weights count repeated rows (n = Σw); reliability weights are
/// inverse variances (the effective n is (Σw)²/Σw²); sampling weights are
/// inverse inclusion probabilities from a survey design.
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum WeightKind {
  Frequency
  Reliability
  Sampling
}

///|
/// Weights with their meaning, without the values they weight: what a
/// regression design or a weighted summary takes next to its data.
#derive(Debug)
pub struct Weights {
  values : FixedArray[Double]
  kind : WeightKind
}

///|
/// Rejects negative and NaN weights.
pub fn Weights::new(
  weights : ArrayView[Double],
  kind~ : WeightKind,
) -> Weights raise DataError {
  for i, w in weights {
    if !(w >= 0.0) {
      raise NegativeWeight(index=i)
    }
  }
  { values: FixedArray::from_iter(weights.iter()), kind, }
}

///|
/// The number of weights.
pub fn Weights::length(self : Weights) -> Int32 {
  self.values.length()
}

///|
/// Values with one weight each.
#derive(Debug)
pub struct Weighted {
  values : FixedArray[Double]
  weights : Weights
}

///|
/// Values paired with one weight each; raises `LengthMismatch` when the
/// counts differ.
pub fn Weighted::new(
  values : ArrayView[Double],
  weights : Weights,
) -> Weighted raise DataError {
  if values.length() != weights.length() {
    raise LengthMismatch(expected=values.length(), got=weights.length())
  }
  { values: FixedArray::from_iter(values.iter()), weights, }
}

///|
/// An r × c table of counts, row-major. Tests of independence, Fisher's
/// exact test and measures of association take this, not raw labels.
#derive(Eq, Debug)
pub struct Table {
  rows : Int32
  cols : Int32
  counts : FixedArray[Int64]
}

///|
/// A table of `rows` × `cols` counts, row-major; the counts must be
/// non-negative and fill the table exactly.
pub fn Table::new(
  rows : Int32,
  cols : Int32,
  counts : ArrayView[Int64],
) -> Table raise DataError {
  if rows <= 0 || cols <= 0 {
    raise InvalidShape(rows~, cols~)
  }
  if counts.length() != rows * cols {
    raise LengthMismatch(expected=rows * cols, got=counts.length())
  }
  for i, c in counts {
    if c < 0 {
      raise NegativeCount(row=i / cols, col=i % cols)
    }
  }
  { rows, cols, counts: FixedArray::from_iter(counts.iter()), }
}

///|
/// Cross-tabulate two label columns of equal length: cell (r, c) counts
/// the positions where the row label is r and the column label is c.
/// Labels are `0 ..< rows` and `0 ..< cols`.
pub fn Table::from_labels(
  row_labels : ArrayView[Int32],
  col_labels : ArrayView[Int32],
  rows~ : Int32,
  cols~ : Int32,
) -> Table raise DataError {
  if rows <= 0 || cols <= 0 {
    raise InvalidShape(rows~, cols~)
  }
  if row_labels.length() != col_labels.length() {
    raise LengthMismatch(expected=row_labels.length(), got=col_labels.length())
  }
  let counts = FixedArray::make(rows * cols, 0L)
  for i in 0..<row_labels.length() {
    let r = row_labels[i]
    let c = col_labels[i]
    if r < 0 || r >= rows {
      raise InvalidLabel(index=i, label=r)
    }
    if c < 0 || c >= cols {
      raise InvalidLabel(index=i, label=c)
    }
    counts[r * cols + c] = counts[r * cols + c] + 1L
  }
  { rows, cols, counts, }
}

///|
/// The count in a cell. Row and column are checked separately, so a column
/// past the end never reads the next row.
pub fn Table::at(self : Table, row : Int32, col : Int32) -> Int64 {
  guard row >= 0 && row < self.rows && col >= 0 && col < self.cols else {
    abort(
      "Table::at: (\{row}, \{col}) is outside a \{self.rows} x \{self.cols} table",
    )
  }
  self.counts[row * self.cols + col]
}

///|
/// The sum of all counts.
pub fn Table::total(self : Table) -> Int64 {
  let mut s = 0L
  for c in self.counts {
    s = s + c
  }
  s
}

///|
/// The margin over columns: one total per row.
pub fn Table::row_totals(self : Table) -> FixedArray[Int64] {
  let out = FixedArray::make(self.rows, 0L)
  for i, c in self.counts {
    out[i / self.cols] = out[i / self.cols] + c
  }
  out
}

///|
/// The margin over rows: one total per column.
pub fn Table::col_totals(self : Table) -> FixedArray[Int64] {
  let out = FixedArray::make(self.cols, 0L)
  for i, c in self.counts {
    out[i % self.cols] = out[i % self.cols] + c
  }
  out
}

///|
/// Apply a missing-value policy to one sample. Returns the cleaned values
/// and the number omitted; `Propagate` keeps NaN in place.
pub fn clean(
  values : ArrayView[Double],
  policy : MissingPolicy,
) -> (FixedArray[Double], Int32) raise DataError {
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
