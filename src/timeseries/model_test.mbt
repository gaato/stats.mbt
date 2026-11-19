///|
test "ARIMA orders and state-space models are checked when built" {
  let o = @timeseries.ArimaOrder::new(p=1, d=1, q=1)
  assert_eq((o.seasonal_p, o.period), (0, 1))
  let seasonal = @timeseries.ArimaOrder::new(
    p=0,
    d=1,
    q=1,
    seasonal_d=1,
    seasonal_q=1,
    period=12,
  )
  assert_eq(seasonal.period, 12)
  let negative = try @timeseries.ArimaOrder::new(p=-1, d=0, q=0) catch {
    @regress.InvalidArgument(name~, ..) => name == "p"
    _ => false
  } noraise {
    _ => false
  }
  assert_true(negative)
  let no_period = try
    @timeseries.ArimaOrder::new(p=0, d=0, q=0, seasonal_p=1)
  catch {
    @regress.InvalidArgument(name~, ..) => name == "period"
    _ => false
  } noraise {
    _ => false
  }
  assert_true(no_period)
  // a local level model: one state, a random walk observed with noise
  let one = @data.DataMatrix::new(1, 1, [1.0])
  let m = @timeseries.StateSpace::new(
    transition=one,
    selection=one,
    observation=one,
    state_noise=one,
    observation_noise=0.5,
  )
  assert_eq(m.observation_noise, 0.5)
  let wide = try
    @timeseries.StateSpace::new(
      transition=one,
      selection=one,
      observation=@data.DataMatrix::new(1, 2, [1.0, 0.0]),
      state_noise=one,
      observation_noise=0.5,
    )
  catch {
    @regress.LengthMismatch(..) => true
    _ => false
  } noraise {
    _ => false
  }
  assert_true(wide)
}
