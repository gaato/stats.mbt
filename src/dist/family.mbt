///|
/// Families: a distribution with its parameters left open.
///
/// MoonBit traits have one `Self` and no associated types, so a trait on a
/// family cannot name the distribution type it produces. A family is
/// therefore a record of functions, generic in the distribution type, like
/// Luna-Flow's container dictionaries. A record also carries what a static
/// trait method could not: fixed hyperparameters (the n of a binomial) and
/// a choice among parameterisations (gamma by rate or by scale), each as its
/// own family value.

///|
/// Where a parameter may live, and how the optimiser sees it. Fitting works
/// on the unconstrained line; each constraint maps there and back.
#non_exhaustive
#derive(Eq, Debug)
pub(all) enum Constraint {
  Real
  Positive
  UnitInterval
  Interval(lo~ : Double, hi~ : Double)
  /// x > bound.
  Above(Double)
  /// x < bound.
  Below(Double)
}

///|
/// Whether x satisfies the constraint.
pub fn Constraint::contains(self : Constraint, x : Double) -> Bool {
  match self {
    Real => !x.is_nan() && !x.is_inf()
    Positive => x > 0.0 && !x.is_inf()
    UnitInterval => x > 0.0 && x < 1.0
    Interval(lo~, hi~) => x > lo && x < hi
    Above(lo) => x > lo && !x.is_inf()
    Below(hi) => x < hi && !x.is_inf()
  }
}

///|
/// Identity, log, logit and scaled logit.
pub fn Constraint::to_unconstrained(self : Constraint, x : Double) -> Double {
  match self {
    Real => x
    Positive => @math.ln(x)
    UnitInterval => @math.ln(x / (1.0 - x))
    Interval(lo~, hi~) => {
      let u = (x - lo) / (hi - lo)
      @math.ln(u / (1.0 - u))
    }
    Above(lo) => @math.ln(x - lo)
    Below(hi) => @math.ln(hi - x)
  }
}

///|
/// The value the constraint maps an unconstrained y to, for optimisers
/// that search all reals.
pub fn Constraint::from_unconstrained(self : Constraint, y : Double) -> Double {
  match self {
    Real => y
    Positive => @math.exp(y)
    UnitInterval => 1.0 / (1.0 + @math.exp(-y))
    Interval(lo~, hi~) => lo + (hi - lo) / (1.0 + @math.exp(-y))
    Above(lo) => lo + @math.exp(y)
    Below(hi) => hi - @math.exp(y)
  }
}

///|
/// A free parameter of a family: its name and where it may lie.
#derive(Eq, Debug)
pub(all) struct ParamSpec {
  name : String
  constraint : Constraint
}

///|
/// The free parameters of a family, a constructor from a parameter vector
/// and the inverse. Fixed hyperparameters are captured by the closures.
/// A record anyone can build, so a distribution defined outside this
/// package gets a family, and with it fitting, the same way.
pub(all) struct Family[D] {
  name : String
  specs : Array[ParamSpec]
  make : (ArrayView[Double]) -> D raise ParamError
  params : (D) -> FixedArray[Double]
}

///|
/// The number of free parameters.
pub fn[D] Family::dim(self : Family[D]) -> Int32 {
  self.specs.length()
}

///|
/// Build a distribution from a point on the unconstrained line, the form
/// an optimiser works in.
pub fn[D] Family::from_unconstrained(
  self : Family[D],
  y : ArrayView[Double],
) -> D raise ParamError {
  let x = FixedArray::make(y.length(), 0.0)
  for i, s in self.specs {
    x[i] = s.constraint.from_unconstrained(y[i])
  }
  (self.make)(x[:])
}

///|
/// The same family with one parameter held at `value` (SciPy's `floc`,
/// `fscale`, `f0`): it leaves the free parameters, so a fit does not move
/// it and its count of parameters drops by one.
/// `Exponential::family().fix("loc", 0.0)` is the one-parameter exponential.
pub fn[D] Family::fix(
  self : Family[D],
  name : String,
  value : Double,
) -> Family[D] raise ParamError {
  let mut at = -1
  for i, s in self.specs {
    if s.name == name {
      at = i
    }
  }
  if at < 0 {
    raise InvalidParam(
      name~,
      value~,
      reason="the family has no free parameter of this name",
    )
  }
  if !self.specs[at].constraint.contains(value) {
    raise InvalidParam(name~, value~, reason="outside the parameter's range")
  }
  let specs = []
  for i, s in self.specs {
    if i != at {
      specs.push(s)
    }
  }
  let make = self.make
  let params = self.params
  {
    name: self.name,
    specs,
    make: p => {
      let full = FixedArray::makei(p.length() + 1, i => {
        if i < at {
          p[i]
        } else if i == at {
          value
        } else {
          p[i - 1]
        }
      })
      make(full[:])
    },
    params: d => {
      let full = params(d)
      FixedArray::makei(full.length() - 1, i => {
        if i < at {
          full[i]
        } else {
          full[i + 1]
        }
      })
    },
  }
}

///|
/// An exponential family: p(x | θ) = h(x) exp(η(θ)·T(x) - A(θ)). Fitting
/// one by maximum likelihood is matching the mean of T, so `fit` can do it
/// for every member without an optimiser, and conjugate priors and GLM
/// links follow from the same pieces.
pub(all) struct ExpFamily[D] {
  family : Family[D]
  /// Length of η and T.
  natural_dim : Int32
  natural : (D) -> FixedArray[Double]
  sufficient : (Double) -> FixedArray[Double]
  log_partition : (D) -> Double
  log_base : (Double) -> Double
  /// The distribution whose E[T] equals the given vector, or an error when
  /// none exists (for example a sample variance of zero).
  from_mean_sufficient : (ArrayView[Double]) -> D raise FitError
}

///|
/// A prior P conjugate to a likelihood D: the posterior after data is again
/// a P. A relation between two types, so a record, not a trait.
pub(all) struct Conjugate[P, D] {
  name : String
  update : (P, ArrayView[Double]) -> P raise ParamError
  /// The likelihood for given parameters drawn from the prior.
  likelihood : (ArrayView[Double]) -> D raise ParamError
}

///|
/// The log-likelihood of independent observations under a density. The
/// sum is compensated, since it often adds thousands of terms of mixed size.
pub fn[D : Density] loglik(d : D, data : ArrayView[Double]) -> Double {
  let mut sum = 0.0
  let mut c = 0.0
  for x in data {
    let y = d.logpdf(x) - c
    let t = sum + y
    c = t - sum - y
    sum = t
  }
  sum
}

///|
/// The log-likelihood of integer data, by a compensated sum of `logpmf`.
pub fn[D : Mass] loglik_int(d : D, data : ArrayView[Int64]) -> Double {
  let mut sum = 0.0
  let mut c = 0.0
  for k in data {
    let y = d.logpmf(k) - c
    let t = sum + y
    c = t - sum - y
    sum = t
  }
  sum
}

///|
/// Gradient of log p(x | θ) with respect to the free parameters of a
/// family, one entry per parameter. Families provide it in closed form; it
/// feeds optimisers and the observed information matrix.
pub(open) trait Score {
  fn score(Self, Double) -> FixedArray[Double]
}

// ---------------------------------------------------------------------------
// Families for the distributions defined so far.
// ---------------------------------------------------------------------------

///|
/// With μ = loc, σ = scale: η = (μ/σ², -1/(2σ²)), T(x) = (x, x²),
/// A = μ²/(2σ²) + ln σ, h(x) = -ln(2π)/2.
pub fn Normal::exp_family() -> ExpFamily[Normal] {
  {
    family: Normal::family(),
    natural_dim: 2,
    natural: d => {
      let v = d.scale * d.scale
      [d.loc / v, -0.5 / v]
    },
    sufficient: x => [x, x * x],
    log_partition: d => {
      d.loc * d.loc / (2.0 * d.scale * d.scale) + @math.ln(d.scale)
    },
    log_base: _ => -0.5 * @math.ln(2.0 * @math.PI),
    from_mean_sufficient: m => {
      let variance = m[1] - m[0] * m[0]
      if !(variance > 0.0) {
        raise OutsideSupport(value=variance)
      }
      Normal::new(loc=m[0], scale=variance.sqrt()) catch {
        _ => raise OutsideSupport(value=variance)
      }
    },
  }
}

///|
/// Normal mean with known variance and a normal prior on the mean: the
/// posterior mean is the precision-weighted average.
pub fn Normal::known_variance_prior(
  sigma : Double,
) -> Conjugate[Normal, Normal] {
  {
    name: "normal-normal (known variance)",
    update: (prior, data) => {
      let n = data.length().to_double()
      let mut s = 0.0
      for x in data {
        s = s + x
      }
      let p0 = 1.0 / (prior.scale * prior.scale)
      let pd = n / (sigma * sigma)
      let post_var = 1.0 / (p0 + pd)
      Normal::new(
        loc=post_var * (p0 * prior.loc + s / (sigma * sigma)),
        scale=post_var.sqrt(),
      )
    },
    likelihood: theta => Normal::new(loc=theta[0], scale=sigma),
  }
}

///|
/// η = logit p, T(x) = x, A = n ln(1 + e^η), h(x) = ln C(n, x).
pub fn Binomial::exp_family(n~ : Int64) -> ExpFamily[Binomial] {
  {
    family: Binomial::family(n~),
    natural_dim: 1,
    natural: d => [@math.ln(d.p / (1.0 - d.p))],
    sufficient: x => [x],
    log_partition: d => -n.to_double() * @math.ln(1.0 - d.p),
    log_base: x => {
      @special.lgamma(n.to_double() + 1.0) -
      @special.lgamma(x + 1.0) -
      @special.lgamma(n.to_double() - x + 1.0)
    },
    from_mean_sufficient: m => {
      let p = m[0] / n.to_double()
      Binomial::new(n~, p~) catch {
        _ => raise OutsideSupport(value=m[0])
      }
    },
  }
}
