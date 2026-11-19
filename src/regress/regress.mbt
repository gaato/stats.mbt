///|
/// Linear and generalised linear models. A model is fitted from a design
/// matrix (no formula language: parsing formulas belongs to a data-frame
/// library), and every fit returns its diagnostics with it. Needs dense
/// linear algebra (QR with column pivoting, Cholesky) for the bodies.

///|
/// Why a model could not be fitted. Shared by the model layers above
/// (survival, time series, multivariate, Bayesian), which fail for the
/// same reasons.
#non_exhaustive
#derive(Debug)
pub(all) suberror ModelError {
  TooFewObservations(needed~ : Int32, got~ : Int32)
  /// Two inputs that must have the same length (or the same number of
  /// rows) do not.
  LengthMismatch(expected~ : Int32, got~ : Int32)
  /// An input value that is NaN.
  ContainsNaN(index~ : Int32)
  /// An argument outside what the procedure accepts (a number of folds
  /// below two, a negative penalty).
  InvalidArgument(name~ : String, reason~ : String)
  /// A matrix that must be invertible or positive definite is not
  /// (collinear regressors, a covariance without full rank).
  Singular
  NotConverged(iterations~ : Int32)
  /// Complete separation in a binary GLM: the maximum likelihood estimate
  /// is at infinity.
  Separation
  BadResponse(index~ : Int32, value~ : Double)
}

///|
/// The regressors, with what the fit should add or hold fixed. Built with
/// `Design::new`, which checks that every per-row input has a value for
/// each row of `x`.
#derive(Debug)
pub struct Design {
  x : @data.DataMatrix
  intercept : Bool
  /// Known additive term on the linear predictor (log exposure in a
  /// Poisson rate model).
  offset : FixedArray[Double]?
  /// One weight per row, with its meaning.
  weights : @data.Weights?
  /// One name per column of `x`; "x0", "x1", ... when not given.
  names : Array[String]
}

///|
pub fn Design::new(
  x : @data.DataMatrix,
  intercept~ : Bool,
  offset? : ArrayView[Double],
  weights? : @data.Weights,
  names? : Array[String],
) -> Design raise ModelError {
  if offset is Some(o) && o.length() != x.rows {
    raise LengthMismatch(expected=x.rows, got=o.length())
  }
  if weights is Some(w) && w.length() != x.rows {
    raise LengthMismatch(expected=x.rows, got=w.length())
  }
  let names = match names {
    Some(ns) => {
      if ns.length() != x.cols {
        raise LengthMismatch(expected=x.cols, got=ns.length())
      }
      ns
    }
    None => Array::makei(x.cols, j => "x\{j}")
  }
  {
    x,
    intercept,
    offset: offset.map(o => FixedArray::from_iter(o.iter())),
    weights,
    names,
  }
}

///|
#derive(Eq, Debug)
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
#derive(Eq, Debug)
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
/// An ordinary, weighted or generalised least-squares fit: everything that
/// rests on the sampling distribution of the estimates.
#derive(Eq, Debug)
pub struct LinearFit {
  /// One name per coefficient, the intercept first when there is one.
  names : Array[String]
  coef : FixedArray[Double]
  se : FixedArray[Double]
  covariance : @data.DataMatrix
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
#non_exhaustive
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
#non_exhaustive
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
/// The response of a GLM. A binomial response is often grouped: successes
/// out of a number of trials per row, which a single column cannot carry.
pub(all) enum GlmResponse {
  /// One value per row: 0/1 for a binary response, a count, a measurement.
  Values(ArrayView[Double])
  /// Successes out of trials per row, for a binomial family.
  Counts(successes~ : ArrayView[Int64], trials~ : ArrayView[Int64])
}

///|
#derive(Eq, Debug)
pub struct GlmFit {
  names : Array[String]
  coef : FixedArray[Double]
  se : FixedArray[Double]
  covariance : @data.DataMatrix
  /// The fitted means, on the response scale.
  fitted : FixedArray[Double]
  deviance : Double
  null_deviance : Double
  df_residual : Int32
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
  y : GlmResponse,
  design : Design,
  family : GlmFamily,
  link : Link,
  max_iterations~ : Int32 = 50,
) -> GlmFit raise ModelError {
  ignore((y, design, family, link, max_iterations))
  ...
}

///|
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum Penalty {
  Ridge(lambda~ : Double)
  Lasso(lambda~ : Double)
  ElasticNet(lambda~ : Double, alpha~ : Double)
}

///|
/// A penalised fit. The estimates are biased by design, so it carries no
/// standard errors, ANOVA table or influence measures: those rest on a
/// sampling distribution the penalty changes.
#derive(Eq, Debug)
pub struct PenalizedFit {
  names : Array[String]
  coef : FixedArray[Double]
  fitted : FixedArray[Double]
  residuals : FixedArray[Double]
  penalty : Penalty
  /// The effective degrees of freedom (the trace of the hat matrix for
  /// ridge, the count of non-zero coefficients for the lasso).
  df : Double
  iterations : Int32
  converged : Bool
}

///|
/// Penalised least squares (coordinate descent for L1).
pub fn penalized(
  y : ArrayView[Double],
  design : Design,
  penalty : Penalty,
) -> PenalizedFit raise ModelError {
  ignore((y, design, penalty))
  ...
}

///|
/// k-fold cross-validated prediction error for a fitting procedure, given
/// as a function from a training set to predictions for a test set.
pub fn cross_validate(
  y : ArrayView[Double],
  design : Design,
  fit_predict : (ArrayView[Double], Design, Design) -> FixedArray[Double] raise ModelError,
  folds~ : Int32,
  rng~ : &@num.RandomSource,
) -> @infer.Estimate raise ModelError {
  ignore((y, design, fit_predict, folds, rng))
  ...
}

///|
/// A nonlinear least-squares fit, with the asymptotic covariance of the
/// parameters from the Jacobian at the solution.
#derive(Eq, Debug)
pub struct NonlinearFit {
  params : FixedArray[Double]
  se : FixedArray[Double]
  covariance : @data.DataMatrix
  fitted : FixedArray[Double]
  residuals : FixedArray[Double]
  /// Residual sum of squares.
  rss : Double
  sigma : Double
  df_residual : Int32
  iterations : Int32
  converged : Bool
}

///|
/// Nonlinear least squares (Levenberg–Marquardt) for a model f(x, θ) with
/// its gradient in θ; the NIST StRD nonlinear problems are the test set.
pub fn nonlinear_least_squares(
  x : ArrayView[Double],
  y : ArrayView[Double],
  model~ : (Double, ArrayView[Double]) -> Double,
  gradient~ : (Double, ArrayView[Double]) -> FixedArray[Double],
  start~ : ArrayView[Double],
) -> NonlinearFit raise ModelError {
  ignore((x, y, model, gradient, start))
  ...
}
