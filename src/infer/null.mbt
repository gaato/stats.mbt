///|
/// The null model: an observed statistic next to its null distribution.

///|
/// Observed statistic and null distribution, before any p-value. Keeping
/// the null distribution lets a caller look at it, plot it, or compute a
/// p-value for another alternative without recomputing.
///
/// Built by `NullModel::new` for a null whose two-sided p-value doubles
/// the smaller tail, and by `NullModel::discrete` or `NullModel::atoms`
/// for one that can also give it by minimum likelihood. That rule needs
/// the null's mass function, which a `&Univariate` does not carry, so the
/// value is worked out when the model is built, from the concrete null.
/// Every model therefore answers every alternative.
pub struct NullModel {
  observed : Double
  null_dist : &@dist.Univariate
  two_sided : TwoSidedRule
  p_kind : PValueKind
  procedure : String
  /// The two-sided p-value by minimum likelihood; present exactly when
  /// `two_sided` is `MinLikelihood`.
  priv min_likelihood : Double?
}

///|
/// A null model whose two-sided p-value is twice the smaller tail: the
/// rule for continuous nulls, and the one randomisation nulls use.
pub fn NullModel::new(
  observed~ : Double,
  null_dist~ : &@dist.Univariate,
  p_kind~ : PValueKind,
  procedure~ : String,
) -> NullModel {
  {
    observed,
    null_dist,
    two_sided: Doubling,
    p_kind,
    procedure,
    min_likelihood: None,
  }
}

///|
/// A null model on the integers with one mode (binomial, hypergeometric,
/// Poisson), with either two-sided rule.
pub fn[
  D : @dist.Univariate + @dist.Mass + @dist.IntCdf + @dist.Support + @dist.Mode,
] NullModel::discrete(
  observed~ : Int64,
  null_dist~ : D,
  two_sided~ : TwoSidedRule,
  p_kind~ : PValueKind,
  procedure~ : String,
) -> NullModel {
  let min_likelihood = match two_sided {
    Doubling => None
    MinLikelihood => Some(min_likelihood_p(null_dist, observed))
  }
  {
    observed: observed.to_double(),
    null_dist,
    two_sided,
    p_kind,
    procedure,
    min_likelihood,
  }
}

///|
/// A null model on finitely many atoms anywhere on the line (an exact
/// randomisation distribution, a rank statistic with ties), with either
/// two-sided rule. The minimum-likelihood p-value adds up every atom no
/// more likely than the observed value; an observed value that is not an
/// atom has probability 0, and so has that p-value.
pub fn NullModel::atoms(
  observed~ : Double,
  null_dist~ : @dist.Atoms,
  two_sided~ : TwoSidedRule,
  p_kind~ : PValueKind,
  procedure~ : String,
) -> NullModel {
  let min_likelihood = match two_sided {
    Doubling => None
    MinLikelihood => {
      let bound = null_dist.mass(observed) * (1.0 + 1.0e-7)
      let mut s = @num.CompensatedSum::zero()
      for q in null_dist.probs() {
        if q <= bound {
          s = s.add(q)
        }
      }
      Some(if s.value() > 1.0 { 1.0 } else { s.value() })
    }
  }
  { observed, null_dist, two_sided, p_kind, procedure, min_likelihood, }
}

///|
/// The p-value from a null model, the same code for theory-based and
/// randomisation tests. Each tail includes the observed value itself
/// (P(T >= t), P(T <= t)), which matters for discrete and permutation nulls,
/// and neither tail is computed as 1 minus the other.
pub fn NullModel::p_value(
  self : NullModel,
  alternative : Alternative,
) -> Double {
  let upper = @dist.Cdf::sf_at_least(self.null_dist, self.observed)
  let lower = @dist.Cdf::cdf(self.null_dist, self.observed)
  match alternative {
    Greater => upper
    Less => lower
    TwoSided =>
      match self.min_likelihood {
        Some(p) => p
        None => {
          let p = 2.0 * (if upper < lower { upper } else { lower })
          if p > 1.0 {
            1.0
          } else {
            p
          }
        }
      }
  }
}

///|
/// What is printed of a null model: everything but the distribution,
/// which is a trait object.
#derive(Debug)
priv struct NullModelSummary {
  observed : Double
  two_sided : TwoSidedRule
  p_kind : PValueKind
  procedure : String
}

///|
pub impl Debug for NullModel with fn to_repr(self) {
  let summary : NullModelSummary = {
    observed: self.observed,
    two_sided: self.two_sided,
    p_kind: self.p_kind,
    procedure: self.procedure,
  }
  Repr(summary)
}

///|
/// An `Int64` from an integer-valued Double, the ends of the range
/// standing for everything beyond them.
fn saturate_int64(x : Double) -> Int64 {
  if x.is_nan() {
    0L
  } else if x >= 9.2233720368547758e18 {
    9223372036854775807L
  } else if x <= -9.2233720368547758e18 {
    -9223372036854775808L
  } else {
    x.to_int64()
  }
}

///|
/// Two-sided p-value of a unimodal discrete null by minimum likelihood: the
/// total probability of outcomes no more likely than the observed k, with
/// the relative slack 1 + 1e-7 that R and SciPy use so that outcomes equal
/// in exact arithmetic are not split by rounding. The support and the mode
/// come from the distribution: the other tail starts at the first point
/// beyond the mode that is no more likely than k.
fn[D : @dist.Mass + @dist.IntCdf + @dist.Support + @dist.Mode] min_likelihood_p(
  d : D,
  k : Int64,
) -> Double {
  let lo = saturate_int64(@dist.Support::support_lower(d))
  let hi = saturate_int64(@dist.Support::support_upper(d))
  let mode = saturate_int64(@dist.Mode::mode(d))
  let bound = @dist.Mass::pmf(d, k) * (1.0 + 1.0e-7)
  let p = if k < mode {
    let mut j = mode
    while j <= hi && @dist.Mass::pmf(d, j) > bound {
      j = j + 1
    }
    @dist.IntCdf::cdf_int(d, k) + @dist.IntCdf::sf_int(d, j - 1)
  } else if k > mode {
    let mut j = mode
    while j >= lo && @dist.Mass::pmf(d, j) > bound {
      j = j - 1
    }
    @dist.IntCdf::cdf_int(d, j) + @dist.IntCdf::sf_int(d, k - 1)
  } else {
    1.0
  }
  if p > 1.0 {
    1.0
  } else {
    p
  }
}
