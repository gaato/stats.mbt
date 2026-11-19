///|
/// Linear and generalised linear models. A model is fitted from a design
/// matrix (no formula language: parsing formulas belongs to a data-frame
/// library), and every fit returns its diagnostics with it. Needs dense
/// linear algebra (QR with column pivoting, Cholesky) for the bodies.

///|
/// The regressors, with what the fit should add or hold fixed.
pub struct Design {
  x : @data.DataMatrix
  intercept : Bool
  /// Known additive term on the linear predictor (log exposure in a
  /// Poisson rate model).
  offset : FixedArray[Double]?
  weights : @data.Weights?
  names : Array[String]
}

///|
pub struct AnovaRow {
  source : String
  df : Double
  sum_sq : Double
  mean_sq : Double
  f : Double?
  p_value : Double?
}

///|
/// Influence and fit diagnostics, all per observation unless noted.
pub struct LinearDiagnostics {
  leverage : FixedArray[Double]
  cooks_distance : FixedArray[Double]
  studentized : FixedArray[Double]
  /// Per regressor.
  vif : FixedArray[Double]
  durbin_watson : Double
  /// Condition number of the scaled design, to flag near-collinearity
  /// (the Longley and Filip cases of NIST StRD).
  condition_number : Double
}

///|
pub struct LinearFit {
  coef : FixedArray[Double]
  se : FixedArray[Double]
  /// k × k, row-major.
  covariance : FixedArray[Double]
  fitted : FixedArray[Double]
  residuals : FixedArray[Double]
  sigma : Double
  r_squared : Double
  adj_r_squared : Double
  df_residual : Int32
  anova : Array[AnovaRow]
  loglik : Double
  aic : Double
  bic : Double
  diagnostics : LinearDiagnostics
  /// The design was rank-deficient; aliased columns have `NaN` coefficients.
  rank : Int32
}

///|
#derive(Debug)
pub(all) suberror ModelError {
  TooFewObservations(needed~ : Int32, got~ : Int32)
  NotConverged(iterations~ : Int32)
  Separation
  BadResponse(index~ : Int32, value~ : Double)
}

///|
/// Ordinary or weighted least squares by QR, never by the normal
/// equations (whose condition number is the square of the design's).
pub fn ols(
  y : ArrayView[Double],
  design : Design,
) -> LinearFit raise ModelError {
  ignore((y, design))
  ...
}

///|
/// Generalised least squares with a known error covariance (n × n).
pub fn gls(
  y : ArrayView[Double],
  design : Design,
  error_covariance : @data.DataMatrix,
) -> LinearFit raise ModelError {
  ignore((y, design, error_covariance))
  ...
}

///|
#derive(Eq, Debug)
pub(all) enum Link {
  Identity
  Log
  Logit
  Probit
  CLogLog
  Inverse
  InverseSquare
  Sqrt
}

///|
/// A GLM family: an exponential family response with its variance
/// function and canonical link.
#derive(Eq, Debug)
pub(all) enum GlmFamily {
  Gaussian
  Binomial
  Poisson
  Gamma
  InverseGaussian
  NegativeBinomial(theta~ : Double)
}

///|
pub struct GlmFit {
  coef : FixedArray[Double]
  se : FixedArray[Double]
  covariance : FixedArray[Double]
  deviance : Double
  null_deviance : Double
  dispersion : Double
  loglik : Double
  aic : Double
  bic : Double
  iterations : Int32
  converged : Bool
  /// Deviance and Pearson residuals.
  deviance_residuals : FixedArray[Double]
  pearson_residuals : FixedArray[Double]
}

///|
/// Iteratively reweighted least squares.
pub fn glm(
  y : ArrayView[Double],
  design : Design,
  family : GlmFamily,
  link : Link,
  max_iterations~ : Int32 = 50,
) -> GlmFit raise ModelError {
  ignore((y, design, family, link, max_iterations))
  ...
}

///|
#derive(Eq, Debug)
pub(all) enum Penalty {
  Ridge(lambda~ : Double)
  Lasso(lambda~ : Double)
  ElasticNet(lambda~ : Double, alpha~ : Double)
}

///|
/// Penalised least squares (coordinate descent for L1).
pub fn penalized(
  y : ArrayView[Double],
  design : Design,
  penalty : Penalty,
) -> LinearFit raise ModelError {
  ignore((y, design, penalty))
  ...
}

///|
/// k-fold cross-validated prediction error for a fitting procedure.
pub fn cross_validate(
  y : ArrayView[Double],
  design : Design,
  fit : (ArrayView[Double], Design) -> LinearFit raise ModelError,
  folds~ : Int32,
  rng~ : &@num.RandomSource,
) -> @infer.Estimate raise ModelError {
  ignore((y, design, fit, folds, rng))
  ...
}

///|
/// Nonlinear least squares (Levenberg–Marquardt) for a model f(x, θ) with
/// its Jacobian; the NIST StRD nonlinear problems are the test set.
pub fn nonlinear_least_squares(
  x : ArrayView[Double],
  y : ArrayView[Double],
  model : (Double, ArrayView[Double]) -> Double,
  jacobian : (Double, ArrayView[Double], FixedArray[Double]) -> Unit,
  start : ArrayView[Double],
) -> LinearFit raise ModelError {
  ignore((x, y, model, jacobian, start))
  ...
}
