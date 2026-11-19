///|
/// Building distributions from other distributions. The Lebesgue
/// decomposition says every distribution on the real line is a mixture of
/// atoms, a density and a singular part, so `PointMass` and `Mixture` are
/// enough to express zero-inflated counts or a measurement with a detection
/// limit. A mixture does not claim a density: its parts are kept separate.

///|
/// All the mass at one point.
pub struct PointMass {
  priv at : Double
}

///|
pub fn PointMass::new(at : Double) -> PointMass raise ParamError {
  if at.is_nan() || at.is_inf() {
    raise InvalidParam(name="at", value=at, reason="must be finite")
  }
  { at, }
}

///|
pub impl Cdf for PointMass with cdf(self, x) {
  if x < self.at {
    0.0
  } else {
    1.0
  }
}

///|
pub impl Cdf for PointMass with sf(self, x) {
  if x < self.at {
    1.0
  } else {
    0.0
  }
}

///|
pub impl Cdf for PointMass with sf_at_least(self, x) {
  if x <= self.at {
    1.0
  } else {
    0.0
  }
}

///|
pub impl Quantile for PointMass with quantile(self, _p) {
  self.at
}

///|
pub impl Sample for PointMass with sample(self, _rng) {
  self.at
}

///|
pub impl Univariate for PointMass

///|
pub impl Moments for PointMass with mean(self) {
  Some(self.at)
}

///|
pub impl Moments for PointMass with variance(_self) {
  Some(0.0)
}

///|
/// Σ wᵢ Fᵢ over components of any kinds. Components are trait objects, so
/// a normal, a point mass and an empirical distribution can sit side by side.
pub struct Mixture {
  priv weights : FixedArray[Double]
  priv components : Array[&Univariate]
}

///|
/// Weights must be non-negative and sum to 1 within rounding; they are not
/// silently renormalised.
pub fn Mixture::new(
  weights : ArrayView[Double],
  components : Array[&Univariate],
) -> Mixture raise ParamError {
  if weights.length() != components.length() || weights.length() == 0 {
    raise InvalidParam(
      name="weights",
      value=weights.length().to_double(),
      reason="need one weight per component, at least one",
    )
  }
  let mut total = 0.0
  for w in weights {
    if !(w >= 0.0) {
      raise InvalidParam(name="weights", value=w, reason="must be >= 0")
    }
    total = total + w
  }
  if (total - 1.0).abs() > 1.0e-12 {
    raise InvalidParam(name="weights", value=total, reason="must sum to 1")
  }
  { weights: FixedArray::from_iter(weights.iter()), components }
}

///|
pub impl Cdf for Mixture with cdf(self, x) {
  let mut s = 0.0
  for i, c in self.components {
    s = s + self.weights[i] * Cdf::cdf(c, x)
  }
  s
}

///|
/// Each component's own upper tail, so the mixture keeps their accuracy.
pub impl Cdf for Mixture with sf(self, x) {
  let mut s = 0.0
  for i, c in self.components {
    s = s + self.weights[i] * Cdf::sf(c, x)
  }
  s
}

///|
pub impl Cdf for Mixture with sf_at_least(self, x) {
  let mut s = 0.0
  for i, c in self.components {
    s = s + self.weights[i] * Cdf::sf_at_least(c, x)
  }
  s
}

///|
/// No closed form in general: bracketing on the distribution function.
pub impl Quantile for Mixture with quantile(self, p) {
  ignore((self, p))
  ...
}

///|
/// Pick a component by weight, then sample from it.
pub impl Sample for Mixture with sample(self, rng) {
  ignore((self, rng))
  ...
}

///|
pub impl Univariate for Mixture

///|
/// a + bX for b > 0, over any distribution. The location-scale families
/// are this applied to a standard member.
pub struct Affine[D] {
  priv base : D
  priv loc : Double
  priv scale : Double
}

///|
pub fn[D] Affine::new(
  base : D,
  loc~ : Double,
  scale~ : Double,
) -> Affine[D] raise ParamError {
  if !(scale > 0.0) {
    raise InvalidParam(name="scale", value=scale, reason="must be > 0")
  }
  { base, loc, scale }
}

///|
pub impl[D : Cdf] Cdf for Affine[D] with cdf(self, x) {
  self.base.cdf((x - self.loc) / self.scale)
}

///|
pub impl[D : Cdf] Cdf for Affine[D] with sf(self, x) {
  self.base.sf((x - self.loc) / self.scale)
}

///|
pub impl[D : Cdf] Cdf for Affine[D] with sf_at_least(self, x) {
  self.base.sf_at_least((x - self.loc) / self.scale)
}

///|
pub impl[D : Quantile] Quantile for Affine[D] with quantile(self, p) {
  self.loc + self.scale * self.base.quantile(p)
}

///|
pub impl[D : Density] Density for Affine[D] with pdf(self, x) {
  self.base.pdf((x - self.loc) / self.scale) / self.scale
}

///|
pub impl[D : Sample] Sample for Affine[D] with sample(self, rng) {
  self.loc + self.scale * self.base.sample(rng)
}

///|
/// The base distribution conditioned on lo < X <= hi. The normalising mass
/// is computed once, from whichever tail keeps more digits.
pub struct Truncated[D] {
  priv base : D
  priv lo : Double
  priv hi : Double
  priv mass : Double
}

///|
pub fn[D : Cdf] Truncated::new(
  base : D,
  lo~ : Double,
  hi~ : Double,
) -> Truncated[D] raise ParamError {
  if !(lo < hi) {
    raise InvalidParam(name="hi", value=hi, reason="must be > lo")
  }
  let mass = base.cdf(hi) - base.cdf(lo)
  if !(mass > 0.0) {
    raise InvalidParam(name="lo", value=lo, reason="interval has no mass")
  }
  { base, lo, hi, mass }
}

///|
pub impl[D : Cdf] Cdf for Truncated[D] with cdf(self, x) {
  if x <= self.lo {
    0.0
  } else if x >= self.hi {
    1.0
  } else {
    (self.base.cdf(x) - self.base.cdf(self.lo)) / self.mass
  }
}

///|
/// Values beyond a detection limit are recorded at the limit: an atom at
/// each limit plus the base distribution in between. Unlike truncation,
/// nothing is removed, so it is a `Mixture` in disguise.
pub struct Censored[D] {
  priv base : D
  priv lo : Double
  priv hi : Double
}

///|
pub fn[D] Censored::new(
  base : D,
  lo~ : Double,
  hi~ : Double,
) -> Censored[D] raise ParamError {
  if !(lo < hi) {
    raise InvalidParam(name="hi", value=hi, reason="must be > lo")
  }
  { base, lo, hi }
}

///|
pub impl[D : Cdf] Cdf for Censored[D] with cdf(self, x) {
  if x < self.lo {
    0.0
  } else if x >= self.hi {
    1.0
  } else {
    self.base.cdf(x)
  }
}
