///|
/// Survival analysis on `SurvivalRecord`s (time, event, entry).

///|
/// A step-function estimate of the survival curve with its pointwise
/// standard errors.
pub struct SurvivalCurve {
  times : FixedArray[Double]
  survival : FixedArray[Double]
  se : FixedArray[Double]
  at_risk : FixedArray[Int32]
  events : FixedArray[Int32]
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
pub fn kaplan_meier(records : ArrayView[@data.SurvivalRecord]) -> SurvivalCurve {
  ignore(records)
  ...
}

///|
pub fn SurvivalCurve::interval(
  self : SurvivalCurve,
  level : Double,
  how : CurveInterval,
) -> (FixedArray[Double], FixedArray[Double]) {
  ignore((self, level, how))
  ...
}

///|
/// Nelson–Aalen cumulative hazard.
pub fn nelson_aalen(records : ArrayView[@data.SurvivalRecord]) -> SurvivalCurve {
  ignore(records)
  ...
}

///|
/// Log-rank test across groups, χ² with groups − 1 df (not a hard-coded
/// exp(-x/2), the bug found on mooncakes).
pub fn log_rank(
  groups : ArrayView[ArrayView[@data.SurvivalRecord]],
) -> @infer.TestResult {
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
pub struct CoxFit {
  coef : FixedArray[Double]
  se : FixedArray[Double]
  hazard_ratio : FixedArray[Double]
  loglik : Double
  converged : Bool
}

///|
/// Cox proportional hazards by partial likelihood.
pub fn cox(
  records : ArrayView[@data.SurvivalRecord],
  covariates : @data.DataMatrix,
  ties : Ties,
) -> CoxFit raise @regress.ModelError {
  ignore((records, covariates, ties))
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
