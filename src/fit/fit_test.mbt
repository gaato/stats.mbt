///|
fn near(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol * b.abs().max(1.0)
}

///|
/// For the normal family the exponential-family fit is the sample mean and
/// the divisor-n standard deviation; the log-likelihood and AIC follow.
test "exponential family fit of a normal" {
  let xs = [2.1, 3.9, 3.0, 5.2, 4.4, 2.8]
  let r = @fit.fit_exp_family(@dist.Normal::exp_family(), xs[:])
  let n = xs.length().to_double()
  let mut s = 0.0
  for x in xs {
    s = s + x
  }
  let mean = s / n
  let mut ss = 0.0
  for x in xs {
    ss = ss + (x - mean) * (x - mean)
  }
  let sigma = (ss / n).sqrt()
  assert_true(near(@dist.Moments::mean(r.dist).unwrap(), mean, 1.0e-14))
  assert_true(near(@dist.Moments::variance(r.dist).unwrap(), sigma * sigma, 1.0e-13))
  // at the MLE, ln L = -n/2 (ln(2π σ²) + 1)
  let expected = -0.5 * n * (@math.ln(2.0 * @math.PI * sigma * sigma) + 1.0)
  assert_true(near(r.loglik, expected, 1.0e-13))
  assert_eq(r.k, 2)
  assert_true(near(r.aic(), 4.0 - 2.0 * expected, 1.0e-13))
  assert_true(r.converged)
}

///|
test "a constant sample has no normal fit" {
  let ok = try @fit.fit_exp_family(@dist.Normal::exp_family(), [1.0, 1.0, 1.0][:]) catch {
    _ => false
  } noraise {
    _ => true
  }
  assert_true(!ok)
}

///|
test "akaike weights" {
  let w = @fit.akaike_weights([
    { family: "a", aic: 10.0, bic: 0.0, k: 1 },
    { family: "b", aic: 12.0, bic: 0.0, k: 2 },
  ])
  let e = @math.exp(-1.0)
  assert_true(near(w[0], 1.0 / (1.0 + e), 1.0e-15))
  assert_true(near(w[1], e / (1.0 + e), 1.0e-15))
}
