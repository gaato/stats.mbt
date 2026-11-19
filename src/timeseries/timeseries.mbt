///|
/// Time series: autocorrelation, spectra, ARIMA, state space.

///|
/// Autocorrelations at lags 0 ..= max_lag (lag 0 is 1).
pub fn acf(ts : @data.TimeSeries, max_lag : Int32) -> FixedArray[Double] {
  ignore((ts, max_lag))
  ...
}

///|
/// Partial autocorrelations by the Durbin–Levinson recursion.
pub fn pacf(ts : @data.TimeSeries, max_lag : Int32) -> FixedArray[Double] {
  ignore((ts, max_lag))
  ...
}

///|
/// Ljung–Box test for autocorrelation up to a lag.
pub fn ljung_box(
  residuals : ArrayView[Double],
  lags : Int32,
  fitted_params~ : Int32,
) -> @infer.TestResult {
  ignore((residuals, lags, fitted_params))
  ...
}

///|
#derive(Eq, Debug)
pub(all) enum Taper {
  None
  Hann
  /// Split cosine bell over the given proportion at each end.
  CosineBell(Double)
}

///|
/// Periodogram by FFT (needs the FFT module).
pub struct Spectrum {
  frequencies : FixedArray[Double]
  power : FixedArray[Double]
}

///|
pub fn periodogram(ts : @data.TimeSeries, taper : Taper) -> Spectrum {
  ignore((ts, taper))
  ...
}

///|
pub struct ArimaOrder {
  p : Int32
  d : Int32
  q : Int32
  seasonal_p : Int32
  seasonal_d : Int32
  seasonal_q : Int32
  period : Int32
}

///|
pub struct ArimaFit {
  order : ArimaOrder
  ar : FixedArray[Double]
  ma : FixedArray[Double]
  sigma2 : Double
  loglik : Double
  aic : Double
  residuals : FixedArray[Double]
}

///|
/// Exact maximum likelihood through the state-space form.
pub fn arima(
  ts : @data.TimeSeries,
  order : ArimaOrder,
) -> ArimaFit raise @regress.ModelError {
  ignore((ts, order))
  ...
}

///|
pub fn ArimaFit::forecast(
  self : ArimaFit,
  steps : Int32,
  level : Double,
) -> (FixedArray[Double], Array[@infer.Interval]) {
  ignore((self, steps, level))
  ...
}

///|
/// A linear Gaussian state-space model x' = T x + R η, y = Z x + ε.
pub struct StateSpace {
  transition : @data.DataMatrix
  selection : @data.DataMatrix
  observation : @data.DataMatrix
  state_noise : @data.DataMatrix
  observation_noise : Double
}

///|
/// Kalman filter: filtered states, their variances and the log-likelihood.
pub fn kalman_filter(
  model : StateSpace,
  y : ArrayView[Double],
) -> (@data.DataMatrix, Double) {
  ignore((model, y))
  ...
}

///|
/// Augmented Dickey–Fuller unit-root test.
pub fn adf_test(ts : @data.TimeSeries, lags : Int32) -> @infer.TestResult {
  ignore((ts, lags))
  ...
}
