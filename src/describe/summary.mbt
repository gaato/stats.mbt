///|
/// The summaries asked for most, on a whole array at once. Each rejects
/// empty input and NaN, like `quantile`; the accumulators are the form to
/// use for streams, pieces and parallel work.

///|
/// The compensated sum of the values and their count, after the checks
/// every summary here shares.
fn checked_sum(data : ArrayView[Double]) -> Double raise @data.DataError {
  if data.length() == 0 {
    raise @data.Empty
  }
  let mut s = @num.CompensatedSum::zero()
  for i, x in data {
    if x.is_nan() {
      raise @data.ContainsNaN(index=i)
    }
    s = s.add(x)
  }
  s.value()
}

///|
/// The sum, compensated (Neumaier), so that values of mixed size do not
/// lose their small terms.
pub fn sum(data : ArrayView[Double]) -> Double raise @data.DataError {
  checked_sum(data)
}

///|
/// The mean, from a compensated sum; raises for empty data or a NaN.
pub fn mean(data : ArrayView[Double]) -> Double raise @data.DataError {
  checked_sum(data) / data.length().to_double()
}

///|
/// The variance with divisor n - ddof, in two passes: deviations from the
/// mean are summed, not squares. The divisor is an argument, not a
/// default: `ddof=1` is the sample variance (R, pandas), `ddof=0` the
/// population variance (NumPy's default).
pub fn variance(
  data : ArrayView[Double],
  ddof~ : Int32,
) -> Double raise @data.DataError {
  let m = mean(data)
  let n = data.length()
  if n - ddof <= 0 {
    raise @data.TooFewValues(needed=ddof + 1, got=n)
  }
  let mut ss = @num.CompensatedSum::zero()
  for x in data {
    ss = ss.add((x - m) * (x - m))
  }
  ss.value() / (n - ddof).to_double()
}

///|
/// The square root of `variance`.
pub fn std_dev(
  data : ArrayView[Double],
  ddof~ : Int32,
) -> Double raise @data.DataError {
  variance(data, ddof~).sqrt()
}

///|
/// The middle value, or the mean of the two middle values for an even
/// count: what quantile definitions 2 and 5 to 9 all give at 1/2.
pub fn median(data : ArrayView[Double]) -> Double raise @data.DataError {
  quantile(data, 0.5, Type7)
}
