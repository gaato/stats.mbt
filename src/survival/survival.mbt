///|
/// Survival analysis on `SurvivalRecord`s (time, event, entry). Failures
/// are `@regress.ModelError`, shared by every model layer.

///|
/// A step-function estimate of the survival curve with its pointwise
/// standard errors (Greenwood).
#derive(Eq, Debug)
pub struct SurvivalCurve {
  times : FixedArray[Double]
  survival : FixedArray[Double]
  se : FixedArray[Double]
  at_risk : FixedArray[Int64]
  events : FixedArray[Int64]
}

///|
/// A step-function estimate of the cumulative hazard with its pointwise
/// standard errors. A different quantity from survival, so a different
/// type: S = exp(-H) only for a continuous distribution.
#derive(Eq, Debug)
pub struct HazardCurve {
  times : FixedArray[Double]
  cumulative_hazard : FixedArray[Double]
  se : FixedArray[Double]
  at_risk : FixedArray[Int64]
  events : FixedArray[Int64]
}

///|
/// How pointwise intervals for the survival curve are formed.
#derive(Eq, Debug)
pub(all) enum CurveInterval {
  Plain
  LogLog
  Log
}

///|
/// Kaplan–Meier with Greenwood standard errors; delayed entry respected.
pub fn kaplan_meier(
  records : ArrayView[@data.SurvivalRecord],
) -> SurvivalCurve raise @regress.ModelError {
  ignore(records)
  ...
}

///|
/// One pointwise interval per time of the curve, each with its level and
/// procedure.
pub fn SurvivalCurve::intervals(
  self : SurvivalCurve,
  level~ : Double,
  how : CurveInterval,
) -> Array[@infer.Interval] raise @regress.ModelError {
  ignore((self, level, how))
  ...
}

///|
/// Nelson–Aalen cumulative hazard.
pub fn nelson_aalen(
  records : ArrayView[@data.SurvivalRecord],
) -> HazardCurve raise @regress.ModelError {
  ignore(records)
  ...
}

///|
/// Log-rank test across groups, χ² with groups − 1 df (not a hard-coded
/// exp(-x/2), the bug found on mooncakes).
pub fn log_rank(
  groups : ArrayView[ArrayView[@data.SurvivalRecord]],
) -> @infer.TestResult raise @regress.ModelError {
  ignore(groups)
  ...
}

///|
#derive(Eq, Debug)
pub(all) enum Ties {
  Breslow
  Efron
}

///|
#derive(Eq, Debug)
pub struct CoxFit {
  names : Array[String]
  coef : FixedArray[Double]
  se : FixedArray[Double]
  covariance : @data.DataMatrix
  hazard_ratio : FixedArray[Double]
  loglik : Double
  /// The cumulative baseline hazard at covariates 0 (Breslow's estimator).
  baseline : HazardCurve
  iterations : Int32
  converged : Bool
}

///|
/// Cox proportional hazards by partial likelihood. `covariates` has one
/// row per record.
pub fn cox(
  records : ArrayView[@data.SurvivalRecord],
  covariates : @data.DataMatrix,
  ties : Ties,
) -> CoxFit raise @regress.ModelError {
  if covariates.rows != records.length() {
    raise LengthMismatch(expected=records.length(), got=covariates.rows)
  }
  ignore(ties)
  ...
}

///|
/// Parametric survival model: a family fitted to censored observations
/// (exact values contribute log f, censored ones log S or log F).
pub fn[D : @dist.Density + @dist.Cdf] fit_censored(
  family : @dist.Family[D],
  observations : ArrayView[@data.Observation],
  start~ : ArrayView[Double],
) -> @fit.MleResult[D] raise @dist.FitError {
  ignore((family, observations, start))
  ...
}
