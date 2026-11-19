///|
/// Input shapes beyond a single array, and the conventions that every
/// library resolves differently and this one makes the caller choose.
/// The shapes are the ones found in R, SciPy, statsmodels and the
/// statistics modules on mooncakes (see DESIGN.md for the survey).
///
/// Counts of observations are `Int64` throughout the library (a published
/// summary or a table can exceed an array length); dimensions, indices and
/// labels are `Int32`.

///|
/// The nine sample quantile definitions of Hyndman and Fan (1996).
/// R's default is `Type7`; most mooncakes modules hard-code it.
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum QuantileMethod {
  Type1
  Type2
  Type3
  Type4
  Type5
  Type6
  Type7
  Type8
  Type9
}

///|
/// How tied values share ranks.
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum RankTies {
  /// Mean of the ranks they span (the statistical default).
  Average
  Min
  Max
  /// Consecutive integers without gaps.
  Dense
  /// Order of appearance; no ties.
  Ordinal
}

///|
/// Paired observations (x, y), for correlation, simple regression and
/// paired comparisons. The two arrays always have the same length.
#derive(Debug)
pub struct Bivariate {
  x : FixedArray[Double]
  y : FixedArray[Double]
}

///|
pub fn Bivariate::new(
  x : ArrayView[Double],
  y : ArrayView[Double],
) -> Bivariate raise DataError {
  if x.length() != y.length() {
    raise LengthMismatch(expected=x.length(), got=y.length())
  }
  { x: FixedArray::from_iter(x.iter()), y: FixedArray::from_iter(y.iter()), }
}

///|
/// The number of pairs.
pub fn Bivariate::length(self : Bivariate) -> Int32 {
  self.x.length()
}

///|
/// Values with a group label each (the long format most data arrives in).
/// Labels are 0 ..< groups.
#derive(Debug)
pub struct Grouped {
  values : FixedArray[Double]
  labels : FixedArray[Int32]
  groups : Int32
}

///|
pub fn Grouped::new(
  values : ArrayView[Double],
  labels : ArrayView[Int32],
  groups~ : Int32,
) -> Grouped raise DataError {
  if groups <= 0 {
    raise InvalidShape(rows=values.length(), cols=groups)
  }
  if values.length() != labels.length() {
    raise LengthMismatch(expected=values.length(), got=labels.length())
  }
  for i, g in labels {
    if g < 0 || g >= groups {
      raise InvalidLabel(index=i, label=g)
    }
  }
  {
    values: FixedArray::from_iter(values.iter()),
    labels: FixedArray::from_iter(labels.iter()),
    groups,
  }
}

///|
/// The values of one group, in input order. Empty for a label no value has.
pub fn Grouped::group(self : Grouped, g : Int32) -> FixedArray[Double] {
  let out = []
  for i, label in self.labels {
    if label == g {
      out.push(self.values[i])
    }
  }
  FixedArray::from_iter(out.iter())
}

///|
/// k related measurements on each of n blocks (subjects), row-major
/// n × k: the input of Friedman's test and repeated-measures ANOVA.
#derive(Debug)
pub struct Blocks {
  blocks : Int32
  treatments : Int32
  values : FixedArray[Double]
}

///|
/// `values` holds block 0's treatments first, then block 1's, and so on.
pub fn Blocks::new(
  blocks : Int32,
  treatments : Int32,
  values : ArrayView[Double],
) -> Blocks raise DataError {
  if blocks <= 0 || treatments <= 0 {
    raise InvalidShape(rows=blocks, cols=treatments)
  }
  if values.length() != blocks * treatments {
    raise LengthMismatch(expected=blocks * treatments, got=values.length())
  }
  { blocks, treatments, values: FixedArray::from_iter(values.iter()), }
}

///|
pub fn Blocks::at(self : Blocks, block : Int32, treatment : Int32) -> Double {
  guard block >= 0 &&
    block < self.blocks &&
    treatment >= 0 &&
    treatment < self.treatments else {
    abort(
      "Blocks::at: (\{block}, \{treatment}) is outside \{self.blocks} x \{self.treatments} blocks",
    )
  }
  self.values[block * self.treatments + treatment]
}

///|
/// A sample known only through its summary, as in published tables.
/// Tests that need nothing more (t tests, z tests) accept it.
#derive(Eq, Debug)
pub(all) struct Summary {
  n : Int64
  mean : Double
  sd : Double
}

///|
/// Successes out of trials, for proportion and binomial tests.
#derive(Eq, Debug)
pub(all) struct Proportion {
  successes : Int64
  trials : Int64
}

///|
/// Counts per category with optional null probabilities, for
/// goodness-of-fit tests; without probabilities the null is uniform.
#derive(Debug)
pub(all) struct CategoryCounts {
  counts : FixedArray[Int64]
  null_probs : FixedArray[Double]?
}

///|
/// One subject in a survival study: the time, whether the event was
/// observed (false means right-censored), and the entry time for
/// left-truncated (delayed-entry) data. Every survival module on mooncakes
/// uses a record like this rather than parallel arrays.
#derive(Eq, Debug)
pub(all) struct SurvivalRecord {
  time : Double
  event : Bool
  entry : Double
}

///|
/// A record observed from time 0 unless `entry` says otherwise, so the
/// common case does not have to spell out the entry time.
pub fn SurvivalRecord::new(
  time : Double,
  event~ : Bool,
  entry~ : Double = 0.0,
) -> SurvivalRecord {
  { time, event, entry, }
}

///|
pub fn SurvivalRecord::to_observation(self : SurvivalRecord) -> Observation {
  if self.event {
    Exact(self.time)
  } else {
    RightCensored(self.time)
  }
}

///|
/// A dense row-major matrix of data (observations × variables). A plain
/// container until a dense linear algebra module provides a matrix type;
/// design matrices, multivariate samples and MCMC draws use it.
#derive(Eq, Debug)
pub struct DataMatrix {
  rows : Int32
  cols : Int32
  values : FixedArray[Double]
}

///|
pub fn DataMatrix::new(
  rows : Int32,
  cols : Int32,
  values : ArrayView[Double],
) -> DataMatrix raise DataError {
  if rows <= 0 || cols <= 0 {
    raise InvalidShape(rows~, cols~)
  }
  if values.length() != rows * cols {
    raise LengthMismatch(expected=rows * cols, got=values.length())
  }
  { rows, cols, values: FixedArray::from_iter(values.iter()), }
}

///|
/// From one array per row; every row must have the length of the first.
pub fn DataMatrix::from_rows(
  rows : ArrayView[ArrayView[Double]],
) -> DataMatrix raise DataError {
  if rows.length() == 0 || rows[0].length() == 0 {
    raise InvalidShape(
      rows=rows.length(),
      cols=if rows.length() == 0 { 0 } else { rows[0].length() },
    )
  }
  let cols = rows[0].length()
  let values = FixedArray::make(rows.length() * cols, 0.0)
  for i, row in rows {
    if row.length() != cols {
      raise LengthMismatch(expected=cols, got=row.length())
    }
    for j, v in row {
      values[i * cols + j] = v
    }
  }
  { rows: rows.length(), cols, values, }
}

///|
/// The entry in row i, column j, each checked against its own bound.
pub fn DataMatrix::at(self : DataMatrix, i : Int32, j : Int32) -> Double {
  guard i >= 0 && i < self.rows && j >= 0 && j < self.cols else {
    abort(
      "DataMatrix::at: (\{i}, \{j}) is outside a \{self.rows} x \{self.cols} matrix",
    )
  }
  self.values[i * self.cols + j]
}

///|
/// Row i (one observation), as a view into the matrix.
pub fn DataMatrix::row(self : DataMatrix, i : Int32) -> ArrayView[Double] {
  guard i >= 0 && i < self.rows else {
    abort("DataMatrix::row: \{i} is outside \{self.rows} rows")
  }
  self.values[i * self.cols:(i + 1) * self.cols]
}

///|
/// Column j (one variable), copied out.
pub fn DataMatrix::column(self : DataMatrix, j : Int32) -> FixedArray[Double] {
  guard j >= 0 && j < self.cols else {
    abort("DataMatrix::column: \{j} is outside \{self.cols} columns")
  }
  FixedArray::makei(self.rows, i => self.values[i * self.cols + j])
}

///|
/// Equally spaced observations. `frequency` is the number per cycle
/// (12 for monthly data with a yearly season), used by seasonal models.
#derive(Debug)
pub(all) struct TimeSeries {
  values : FixedArray[Double]
  start : Double
  frequency : Int32
}

///|
/// Angles in radians, for directional statistics.
#derive(Debug)
pub(all) struct Circular {
  angles : FixedArray[Double]
}
