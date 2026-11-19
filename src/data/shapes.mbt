///|
/// Input shapes beyond a single array, and the conventions that every
/// library resolves differently and this one makes the caller choose.
/// The shapes are the ones found in R, SciPy, statsmodels and the
/// statistics modules on mooncakes (see DESIGN.md for the survey).

///|
/// The nine sample quantile definitions of Hyndman and Fan (1996).
/// R's default is `Type7`; most mooncakes modules hard-code it.
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
/// Paired observations (x, y), for correlation and simple regression.
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
/// Values with a group label each (the long format most data arrives in).
/// Labels are 0 ..< groups.
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
  if values.length() != labels.length() {
    raise LengthMismatch(expected=values.length(), got=labels.length())
  }
  ignore(groups)
  ...
}

///|
/// The values of one group, in input order.
pub fn Grouped::group(self : Grouped, g : Int32) -> FixedArray[Double] {
  ignore((self, g))
  ...
}

///|
/// k related measurements on each of n blocks (subjects), row-major
/// n × k: the input of Friedman's test and repeated-measures ANOVA.
pub struct Blocks {
  blocks : Int32
  treatments : Int32
  values : FixedArray[Double]
}

///|
/// A sample known only through its summary, as in published tables.
/// Tests that need nothing more (t tests, z tests) accept it.
pub(all) struct Summary {
  n : Int32
  mean : Double
  sd : Double
}

///|
/// Successes out of trials, for proportion and binomial tests.
pub(all) struct Proportion {
  successes : Int64
  trials : Int64
}

///|
/// Counts per category with optional null probabilities, for
/// goodness-of-fit tests; without probabilities the null is uniform.
pub(all) struct CategoryCounts {
  counts : FixedArray[Int64]
  null_probs : FixedArray[Double]?
}

///|
/// One subject in a survival study: the time, whether the event was
/// observed (false means right-censored), and the entry time for
/// left-truncated (delayed-entry) data. Every survival module on mooncakes
/// uses a record like this rather than parallel arrays.
pub(all) struct SurvivalRecord {
  time : Double
  event : Bool
  entry : Double
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
  if values.length() != rows * cols {
    raise LengthMismatch(expected=rows * cols, got=values.length())
  }
  { rows, cols, values: FixedArray::from_iter(values.iter()), }
}

///|
pub fn DataMatrix::at(self : DataMatrix, i : Int32, j : Int32) -> Double {
  self.values[i * self.cols + j]
}

///|
/// Equally spaced observations. `frequency` is the number per cycle
/// (12 for monthly data with a yearly season), used by seasonal models.
pub(all) struct TimeSeries {
  values : FixedArray[Double]
  start : Double
  frequency : Int32
}

///|
/// Angles in radians, for directional statistics.
pub(all) struct Circular {
  angles : FixedArray[Double]
}
