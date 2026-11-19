///|
/// Time series: autocorrelation, spectra, ARIMA, state space. Functions
/// that only read the values take an array; those that need the sampling
/// frequency (seasonal models, spectra in cycles per unit time) take a
/// `TimeSeries`. Failures are `@regress.ModelError`.

///|
/// Autocorrelations at lags 0 ..= max_lag (lag 0 is 1).
pub fn acf(
  values : ArrayView[Double],
  max_lag~ : Int32,
) -> FixedArray[Double] raise @regress.ModelError {
  ignore((values, max_lag))
  ...
}

///|
/// Partial autocorrelations at lags 1 ..= max_lag, by the Durbin–Levinson
/// recursion.
pub fn pacf(
  values : ArrayView[Double],
  max_lag~ : Int32,
) -> FixedArray[Double] raise @regress.ModelError {
  ignore((values, max_lag))
  ...
}

///|
/// Ljung–Box test for autocorrelation up to a lag.
pub fn ljung_box(
  residuals : ArrayView[Double],
  lags~ : Int32,
  fitted_params~ : Int32,
) -> @infer.TestResult raise @regress.ModelError {
  ignore((residuals, lags, fitted_params))
  ...
}

///|
/// The window applied before a periodogram.
#derive(Eq, Debug)
pub(all) enum Taper {
  /// No taper: every value has weight 1.
  Rectangular
  Hann
  /// Split cosine bell over the given proportion at each end.
  CosineBell(Double)
}

///|
/// Periodogram by FFT (needs the FFT module).
#derive(Eq, Debug)
pub struct Spectrum {
  frequencies : FixedArray[Double]
  power : FixedArray[Double]
}

///|
pub fn periodogram(
  ts : @data.TimeSeries,
  taper : Taper,
) -> Spectrum raise @regress.ModelError {
  ignore((ts, taper))
  ...
}

///|
/// The orders of a seasonal ARIMA(p, d, q)(P, D, Q)ₛ model.
#derive(Eq, Debug)
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
/// The seasonal orders default to 0, which leaves a plain ARIMA(p, d, q).
/// Every order is non-negative, and a seasonal part needs a period of at
/// least 2.
pub fn ArimaOrder::new(
  p~ : Int32,
  d~ : Int32,
  q~ : Int32,
  seasonal_p~ : Int32 = 0,
  seasonal_d~ : Int32 = 0,
  seasonal_q~ : Int32 = 0,
  period~ : Int32 = 1,
) -> ArimaOrder raise @regress.ModelError {
  for
    pair in [
      ("p", p),
      ("d", d),
      ("q", q),
      ("seasonal_p", seasonal_p),
      ("seasonal_d", seasonal_d),
      ("seasonal_q", seasonal_q),
    ] {
    if pair.1 < 0 {
      raise InvalidArgument(name=pair.0, reason="must be >= 0")
    }
  }
  let seasonal = seasonal_p + seasonal_d + seasonal_q > 0
  if period < 1 || (seasonal && period < 2) {
    raise InvalidArgument(
      name="period",
      reason="must be >= 1, and >= 2 with a seasonal part",
    )
  }
  { p, d, q, seasonal_p, seasonal_d, seasonal_q, period, }
}

///|
#derive(Eq, Debug)
pub struct ArimaFit {
  order : ArimaOrder
  ar : FixedArray[Double]
  ma : FixedArray[Double]
  seasonal_ar : FixedArray[Double]
  seasonal_ma : FixedArray[Double]
  /// The mean (or drift, after differencing); `None` when not estimated.
  intercept : Double?
  /// Standard errors of the coefficients, in the order ar, ma, seasonal ar,
  /// seasonal ma, intercept.
  se : FixedArray[Double]
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
/// Point forecasts with prediction intervals.
#derive(Eq, Debug)
pub struct Forecast {
  mean : FixedArray[Double]
  intervals : Array[@infer.Interval]
}

///|
pub fn ArimaFit::forecast(
  self : ArimaFit,
  steps~ : Int32,
  level~ : Double,
) -> Forecast raise @regress.ModelError {
  ignore((self, steps, level))
  ...
}

///|
/// A linear Gaussian state-space model x' = T x + R η, y = Z x + ε.
#derive(Eq, Debug)
pub struct StateSpace {
  transition : @data.DataMatrix
  selection : @data.DataMatrix
  observation : @data.DataMatrix
  state_noise : @data.DataMatrix
  observation_noise : Double
}

///|
/// Checks that the matrices fit together: T is m × m, R is m × r, Z is
/// 1 × m and the state noise covariance r × r.
pub fn StateSpace::new(
  transition~ : @data.DataMatrix,
  selection~ : @data.DataMatrix,
  observation~ : @data.DataMatrix,
  state_noise~ : @data.DataMatrix,
  observation_noise~ : Double,
) -> StateSpace raise @regress.ModelError {
  let m = transition.rows
  if transition.cols != m {
    raise LengthMismatch(expected=m, got=transition.cols)
  }
  if selection.rows != m {
    raise LengthMismatch(expected=m, got=selection.rows)
  }
  if observation.rows != 1 || observation.cols != m {
    raise LengthMismatch(expected=m, got=observation.cols)
  }
  let r = selection.cols
  if state_noise.rows != r || state_noise.cols != r {
    raise LengthMismatch(expected=r, got=state_noise.rows)
  }
  if !(observation_noise >= 0.0) {
    raise InvalidArgument(name="observation_noise", reason="must be >= 0")
  }
  { transition, selection, observation, state_noise, observation_noise, }
}

///|
/// What the Kalman filter returns: the filtered state means (one row per
/// time), the variance of each prediction error, and the log-likelihood.
#derive(Eq, Debug)
pub struct KalmanResult {
  states : @data.DataMatrix
  prediction_variance : FixedArray[Double]
  loglik : Double
}

///|
pub fn kalman_filter(
  model : StateSpace,
  y : ArrayView[Double],
) -> KalmanResult raise @regress.ModelError {
  ignore((model, y))
  ...
}

///|
/// Augmented Dickey–Fuller unit-root test.
pub fn adf_test(
  values : ArrayView[Double],
  lags~ : Int32,
) -> @infer.TestResult raise @regress.ModelError {
  ignore((values, lags))
  ...
}
