///|
/// Building distributions from other distributions. The Lebesgue
/// decomposition says every distribution on the real line is a mixture of
/// atoms, a density and a singular part, so `PointMass` and `Mixture` are
/// enough to express zero-inflated counts or a measurement with a detection
/// limit. A mixture does not claim a density: its parts are kept separate.
///
/// Every wrapper forwards each tail to the same tail of its base. None
/// relies on a trait default, so wrapping a distribution never costs it
/// the digits its own formulas keep.

///|
/// All the mass at one point.
#derive(Eq, Debug)
pub struct PointMass {
  at : Double
}

///|
/// All mass at one finite point.
pub fn PointMass::new(at : Double) -> PointMass raise ParamError {
  if at.is_nan() || at.is_inf() {
    raise InvalidParam(name="at", value=at, reason="must be finite")
  }
  { at, }
}

///|
/// P(X ≤ x).
pub impl Cdf for PointMass with fn cdf(self, x) {
  if x < self.at {
    0.0
  } else {
    1.0
  }
}

///|
/// P(X > x).
pub impl Cdf for PointMass with fn sf(self, x) {
  if x < self.at {
    1.0
  } else {
    0.0
  }
}

///|
/// P(X ≥ x).
pub impl Cdf for PointMass with fn sf_at_least(self, x) {
  if x <= self.at {
    1.0
  } else {
    0.0
  }
}

///|
/// The smallest x with P(X ≤ x) ≥ p.
pub impl Quantile for PointMass with fn quantile(self, _p) {
  self.at
}

///|
/// The smallest x with P(X > x) ≤ q.
pub impl Quantile for PointMass with fn isf(self, _q) {
  self.at
}

///|
/// One draw from `rng`.
pub impl Sample for PointMass with fn sample(self, _rng) {
  self.at
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl Univariate for PointMass

///|
/// The mean, or `None` when it does not exist.
pub impl Moments for PointMass with fn mean(self) {
  Some(self.at)
}

///|
/// The variance, or `None` when it does not exist.
pub impl Moments for PointMass with fn variance(_self) {
  Some(0.0)
}

///|
/// The lower end of the support.
pub impl Support for PointMass with fn support_lower(self) {
  self.at
}

///|
/// The upper end of the support.
pub impl Support for PointMass with fn support_upper(self) {
  self.at
}

///|
/// The most likely value.
pub impl Mode for PointMass with fn mode(self) {
  self.at
}

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend PointMass with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend PointMass with Quantile::{quantile, isf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend PointMass with Sample::{sample, sample_into}

///|
/// Dot syntax for `Moments`: `d.mean`, `d.variance`.
pub extend PointMass with Moments::{mean, variance}

///|
/// Dot syntax for `Support`: `d.support_lower`, `d.support_upper`.
pub extend PointMass with Support::{support_lower, support_upper}

///|
/// Dot syntax for `Mode`: `d.mode`.
pub extend PointMass with Mode::{mode}

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
  { weights: FixedArray::from_iter(weights.iter()), components, }
}

///|
/// The mixing weights, in the order the components were given.
pub fn Mixture::weights(self : Mixture) -> ArrayView[Double] {
  self.weights[:]
}

///|
/// P(X ≤ x).
pub impl Cdf for Mixture with fn cdf(self, x) {
  let mut s = 0.0
  for i, c in self.components {
    s = s + self.weights[i] * Cdf::cdf(c, x)
  }
  s
}

///|
/// Each component's own upper tail, so the mixture keeps their accuracy.
pub impl Cdf for Mixture with fn sf(self, x) {
  let mut s = 0.0
  for i, c in self.components {
    s = s + self.weights[i] * Cdf::sf(c, x)
  }
  s
}

///|
/// P(X ≥ x).
pub impl Cdf for Mixture with fn sf_at_least(self, x) {
  let mut s = 0.0
  for i, c in self.components {
    s = s + self.weights[i] * Cdf::sf_at_least(c, x)
  }
  s
}

///|
/// The smallest x in [lo, hi] where `reached` holds, for a condition that
/// is false below some point and true from it on, and true at `hi`.
/// Bisection down to neighbouring Doubles.
fn smallest_where(
  reached : (Double) -> Bool,
  lo : Double,
  hi : Double,
) -> Double {
  if reached(lo) {
    return lo
  }
  let mut a = lo
  let mut b = hi
  while true {
    let m = a + 0.5 * (b - a)
    if !(m > a && m < b) {
      break
    }
    if reached(m) {
      b = m
    } else {
      a = m
    }
  }
  b
}

///|
/// The range of one quantile over the components that carry weight. The
/// mixture's own quantile lies inside it: below the smallest component
/// quantile every component is short of the probability, and at the
/// largest every component has reached it.
fn Mixture::bracket(
  self : Mixture,
  of : (&Univariate) -> Double,
) -> (Double, Double) {
  let mut lo = 1.0 / 0.0
  let mut hi = -1.0 / 0.0
  for i, c in self.components {
    if self.weights[i] > 0.0 {
      let v = of(c)
      if v < lo {
        lo = v
      }
      if v > hi {
        hi = v
      }
    }
  }
  (lo, hi)
}

///|
/// No closed form in general: bisection on the distribution function,
/// between the component quantiles.
pub impl Quantile for Mixture with fn quantile(self, p) {
  if !(p >= 0.0 && p <= 1.0) {
    return 0.0 / 0.0
  }
  let (lo, hi) = self.bracket(c => Quantile::quantile(c, p))
  smallest_where(x => Cdf::cdf(self, x) >= p, lo, hi)
}

///|
/// The same on the upper tail, so a small q is never turned into 1 - q.
pub impl Quantile for Mixture with fn isf(self, q) {
  if !(q >= 0.0 && q <= 1.0) {
    return 0.0 / 0.0
  }
  let (lo, hi) = self.bracket(c => Quantile::isf(c, q))
  smallest_where(x => Cdf::sf(self, x) <= q, lo, hi)
}

///|
/// Pick a component by weight with one uniform, then sample from it.
pub impl Sample for Mixture with fn sample(self, rng) {
  let u = @num.uniform(rng)
  let mut acc = 0.0
  let mut chosen = 0
  for i, w in self.weights {
    if w > 0.0 {
      chosen = i
      acc = acc + w
      if u < acc {
        break
      }
    }
  }
  Sample::sample(self.components[chosen], rng)
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl Univariate for Mixture

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend Mixture with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend Mixture with Quantile::{quantile, isf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend Mixture with Sample::{sample, sample_into}

///|
/// loc + scale·X for scale > 0, over any distribution. The catalogue's
/// continuous families carry their own `loc` and `scale`; this is for the
/// rest (an empirical distribution, a mixture, a caller's own type).
pub struct Affine[D] {
  base : D
  loc : Double
  scale : Double
}

///|
/// `base` moved to `loc` and stretched by `scale` (positive): X = loc +
/// scale·Y, with every tail forwarded.
pub fn[D] Affine::new(
  base : D,
  loc~ : Double = 0.0,
  scale~ : Double = 1.0,
) -> Affine[D] raise ParamError {
  if loc.is_nan() || loc.is_inf() {
    raise InvalidParam(name="loc", value=loc, reason="must be finite")
  }
  if !(scale > 0.0) || scale.is_inf() {
    raise InvalidParam(name="scale", value=scale, reason="must be in (0, inf)")
  }
  { base, loc, scale, }
}

///|
/// P(X ≤ x).
pub impl[D : Cdf] Cdf for Affine[D] with fn cdf(self, x) {
  self.base.cdf((x - self.loc) / self.scale)
}

///|
/// P(X > x).
pub impl[D : Cdf] Cdf for Affine[D] with fn sf(self, x) {
  self.base.sf((x - self.loc) / self.scale)
}

///|
/// ln P(X ≤ x).
pub impl[D : Cdf] Cdf for Affine[D] with fn logcdf(self, x) {
  self.base.logcdf((x - self.loc) / self.scale)
}

///|
/// ln P(X > x).
pub impl[D : Cdf] Cdf for Affine[D] with fn logsf(self, x) {
  self.base.logsf((x - self.loc) / self.scale)
}

///|
/// P(X ≥ x).
pub impl[D : Cdf] Cdf for Affine[D] with fn sf_at_least(self, x) {
  self.base.sf_at_least((x - self.loc) / self.scale)
}

///|
/// The smallest x with P(X ≤ x) ≥ p.
pub impl[D : Quantile] Quantile for Affine[D] with fn quantile(self, p) {
  self.loc + self.scale * self.base.quantile(p)
}

///|
/// The smallest x with P(X > x) ≤ q.
pub impl[D : Quantile] Quantile for Affine[D] with fn isf(self, q) {
  self.loc + self.scale * self.base.isf(q)
}

///|
/// The density at x.
pub impl[D : Density] Density for Affine[D] with fn pdf(self, x) {
  self.base.pdf((x - self.loc) / self.scale) / self.scale
}

///|
/// The log density at x.
pub impl[D : Density] Density for Affine[D] with fn logpdf(self, x) {
  self.base.logpdf((x - self.loc) / self.scale) - @math.ln(self.scale)
}

///|
/// One draw from `rng`.
pub impl[D : Sample] Sample for Affine[D] with fn sample(self, rng) {
  self.loc + self.scale * self.base.sample(rng)
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl[D : Univariate] Univariate for Affine[D]

///|
/// The mean, or `None` when it does not exist.
pub impl[D : Moments] Moments for Affine[D] with fn mean(self) {
  match self.base.mean() {
    Some(m) => Some(self.loc + self.scale * m)
    None => None
  }
}

///|
/// The variance, or `None` when it does not exist.
pub impl[D : Moments] Moments for Affine[D] with fn variance(self) {
  match self.base.variance() {
    Some(v) => Some(self.scale * self.scale * v)
    None => None
  }
}

///|
/// The lower end of the support.
pub impl[D : Support] Support for Affine[D] with fn support_lower(self) {
  self.loc + self.scale * self.base.support_lower()
}

///|
/// The upper end of the support.
pub impl[D : Support] Support for Affine[D] with fn support_upper(self) {
  self.loc + self.scale * self.base.support_upper()
}

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend Affine with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend Affine with Quantile::{quantile, isf}

///|
/// Dot syntax for `Density`: `d.pdf`, `d.logpdf`.
pub extend Affine with Density::{pdf, logpdf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend Affine with Sample::{sample, sample_into}

///|
/// Dot syntax for `Moments`: `d.mean`, `d.variance`.
pub extend Affine with Moments::{mean, variance}

///|
/// Dot syntax for `Support`: `d.support_lower`, `d.support_upper`.
pub extend Affine with Support::{support_lower, support_upper}

///|
/// P(a < X <= b), as a difference of the tail that keeps its digits: the
/// upper one when a is above the median, the lower one otherwise. A
/// difference of distribution functions next to 1 would be all rounding.
fn[D : Cdf] mass_between(base : D, a : Double, b : Double) -> Double {
  if base.cdf(a) > 0.5 {
    base.sf(a) - base.sf(b)
  } else {
    base.cdf(b) - base.cdf(a)
  }
}

///|
/// The base distribution conditioned on lo < X <= hi; either bound may be
/// left out. The normalising mass and the tails beyond the bounds are
/// computed once.
pub struct Truncated[D] {
  base : D
  lo : Double
  hi : Double
  priv mass : Double
  /// P(X <= lo) and P(X > hi) of the base.
  priv below : Double
  priv above : Double
}

///|
/// `base` conditioned on `lo` ≤ X ≤ `hi`; raises when the interval has no
/// mass.
pub fn[D : Cdf] Truncated::new(
  base : D,
  lo~ : Double = -1.0 / 0.0,
  hi~ : Double = 1.0 / 0.0,
) -> Truncated[D] raise ParamError {
  if !(lo < hi) {
    raise InvalidParam(name="hi", value=hi, reason="must be > lo")
  }
  let mass = mass_between(base, lo, hi)
  if !(mass > 0.0) {
    raise InvalidParam(name="lo", value=lo, reason="interval has no mass")
  }
  { base, lo, hi, mass, below: base.cdf(lo), above: base.sf(hi), }
}

///|
/// P(X ≤ x).
pub impl[D : Cdf] Cdf for Truncated[D] with fn cdf(self, x) {
  if x <= self.lo {
    0.0
  } else if x >= self.hi {
    1.0
  } else {
    mass_between(self.base, self.lo, x) / self.mass
  }
}

///|
/// P(X > x).
pub impl[D : Cdf] Cdf for Truncated[D] with fn sf(self, x) {
  if x <= self.lo {
    1.0
  } else if x >= self.hi {
    0.0
  } else {
    mass_between(self.base, x, self.hi) / self.mass
  }
}

///|
/// The upper tail plus the base's own atom at x, if it has one.
pub impl[D : Cdf] Cdf for Truncated[D] with fn sf_at_least(self, x) {
  if x <= self.lo {
    1.0
  } else if x > self.hi {
    0.0
  } else {
    let atom = self.base.sf_at_least(x) - self.base.sf(x)
    (mass_between(self.base, x, self.hi) + atom) / self.mass
  }
}

///|
fn[D] Truncated::clamp(self : Truncated[D], x : Double) -> Double {
  if x < self.lo {
    self.lo
  } else if x > self.hi {
    self.hi
  } else {
    x
  }
}

///|
/// Inverts the base on whichever tail holds the target probability: the
/// lower target `below + p·mass` when it is at most a half, the upper one
/// `above + (1 - p)·mass` otherwise.
pub impl[D : Quantile] Quantile for Truncated[D] with fn quantile(self, p) {
  if !(p >= 0.0 && p <= 1.0) {
    return 0.0 / 0.0
  }
  let lower = self.below + p * self.mass
  let x = if lower <= 0.5 {
    self.base.quantile(lower)
  } else {
    self.base.isf(self.above + (1.0 - p) * self.mass)
  }
  self.clamp(x)
}

///|
/// The smallest x with P(X > x) ≤ q.
pub impl[D : Quantile] Quantile for Truncated[D] with fn isf(self, q) {
  if !(q >= 0.0 && q <= 1.0) {
    return 0.0 / 0.0
  }
  let upper = self.above + q * self.mass
  let x = if upper <= 0.5 {
    self.base.isf(upper)
  } else {
    self.base.quantile(self.below + (1.0 - q) * self.mass)
  }
  self.clamp(x)
}

///|
/// The density at x.
pub impl[D : Density] Density for Truncated[D] with fn pdf(self, x) {
  if x <= self.lo || x > self.hi {
    0.0
  } else {
    self.base.pdf(x) / self.mass
  }
}

///|
/// The log density at x.
pub impl[D : Density] Density for Truncated[D] with fn logpdf(self, x) {
  if x <= self.lo || x > self.hi {
    -1.0 / 0.0
  } else {
    self.base.logpdf(x) - @math.ln(self.mass)
  }
}

///|
/// Inversion of the quantile function.
pub impl[D : Quantile] Sample for Truncated[D] with fn sample(self, rng) {
  Quantile::quantile(self, @num.uniform_open(rng))
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl[D : Univariate] Univariate for Truncated[D]

///|
/// The lower end of the support.
pub impl[D : Support] Support for Truncated[D] with fn support_lower(self) {
  let a = self.base.support_lower()
  if a > self.lo {
    a
  } else {
    self.lo
  }
}

///|
/// The upper end of the support.
pub impl[D : Support] Support for Truncated[D] with fn support_upper(self) {
  let b = self.base.support_upper()
  if b < self.hi {
    b
  } else {
    self.hi
  }
}

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend Truncated with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend Truncated with Quantile::{quantile, isf}

///|
/// Dot syntax for `Density`: `d.pdf`, `d.logpdf`.
pub extend Truncated with Density::{pdf, logpdf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend Truncated with Sample::{sample, sample_into}

///|
/// Dot syntax for `Support`: `d.support_lower`, `d.support_upper`.
pub extend Truncated with Support::{support_lower, support_upper}

///|
/// Values beyond a detection limit are recorded at the limit: an atom at
/// each limit plus the base distribution in between. Unlike truncation,
/// nothing is removed, so it is a `Mixture` in disguise. Either limit may
/// be left out.
pub struct Censored[D] {
  base : D
  lo : Double
  hi : Double
}

///|
/// `base` clipped to [`lo`, `hi`]: the mass outside moves to the ends.
pub fn[D] Censored::new(
  base : D,
  lo~ : Double = -1.0 / 0.0,
  hi~ : Double = 1.0 / 0.0,
) -> Censored[D] raise ParamError {
  if !(lo < hi) {
    raise InvalidParam(name="hi", value=hi, reason="must be > lo")
  }
  { base, lo, hi, }
}

///|
/// P(X ≤ x).
pub impl[D : Cdf] Cdf for Censored[D] with fn cdf(self, x) {
  if x < self.lo {
    0.0
  } else if x >= self.hi {
    1.0
  } else {
    self.base.cdf(x)
  }
}

///|
/// P(X > x).
pub impl[D : Cdf] Cdf for Censored[D] with fn sf(self, x) {
  if x < self.lo {
    1.0
  } else if x >= self.hi {
    0.0
  } else {
    self.base.sf(x)
  }
}

///|
/// The atom at each limit counts: P(X >= lo) = 1, and P(X >= hi) is the
/// base's mass at or beyond the upper limit.
pub impl[D : Cdf] Cdf for Censored[D] with fn sf_at_least(self, x) {
  if x <= self.lo {
    1.0
  } else if x > self.hi {
    0.0
  } else {
    self.base.sf_at_least(x)
  }
}

///|
fn[D] Censored::clamp(self : Censored[D], x : Double) -> Double {
  if x < self.lo {
    self.lo
  } else if x > self.hi {
    self.hi
  } else {
    x
  }
}

///|
/// The smallest x with P(X ≤ x) ≥ p.
pub impl[D : Quantile] Quantile for Censored[D] with fn quantile(self, p) {
  self.clamp(self.base.quantile(p))
}

///|
/// The smallest x with P(X > x) ≤ q.
pub impl[D : Quantile] Quantile for Censored[D] with fn isf(self, q) {
  self.clamp(self.base.isf(q))
}

///|
/// One draw from `rng`.
pub impl[D : Sample] Sample for Censored[D] with fn sample(self, rng) {
  self.clamp(self.base.sample(rng))
}

///|
/// Usable wherever a `&Univariate` is taken.
pub impl[D : Univariate] Univariate for Censored[D]

///|
/// Dot syntax for `Cdf`: `d.cdf`, `d.sf`, `d.logcdf`, `d.logsf`, `d.sf_at_least`.
pub extend Censored with Cdf::{cdf, sf, logcdf, logsf, sf_at_least}

///|
/// Dot syntax for `Quantile`: `d.quantile`, `d.isf`.
pub extend Censored with Quantile::{quantile, isf}

///|
/// Dot syntax for `Sample`: `d.sample`, `d.sample_into`.
pub extend Censored with Sample::{sample, sample_into}
